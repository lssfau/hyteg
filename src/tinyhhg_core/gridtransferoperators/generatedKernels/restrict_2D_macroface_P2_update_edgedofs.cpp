
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void restrict_2D_macroface_P2_update_edgedofs_level_2(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 11];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 20*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 10];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_3(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 19];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 36*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 18];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_4(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 35];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 68*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 34];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_5(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 67];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 132*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 66];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_6(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 131];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 260*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 130];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_7(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 259];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 516*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 258];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_8(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 515];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 1028*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 514];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_9(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1027];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 2052*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1026];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_10(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2051];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 4100*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2050];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_11(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4099];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 8196*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4098];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_12(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8195];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 16388*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8194];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_13(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16387];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 32772*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16386];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_14(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32771];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + 65540*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32770];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}

static void restrict_2D_macroface_P2_update_edgedofs_level_any(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_35 = 1 / (num_neighbor_faces_edge0);
   const double xi_38 = 1 / (num_neighbor_faces_edge2);
   const double xi_159 = 1 / (num_neighbor_faces_edge0);
   const double xi_66 = 1 / (num_neighbor_faces_edge0);
   const double xi_69 = 1 / (num_neighbor_faces_edge1);
   const double xi_221 = 1 / (num_neighbor_faces_edge2);
   const double xi_190 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_100 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_40 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_56 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_49 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_46 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_47 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_57 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_50 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_51 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_52 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_53 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_54 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_55 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_43 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_44 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_45 = xi_35*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_58 = 1.0*xi_38*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_59 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_60 = xi_38*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_40 + xi_43 + xi_44 + xi_45 + xi_49 + xi_56;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_46 + xi_47 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_57;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)) - 1; ctr_1 += 1)
      {
         const double xi_161 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_180 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_170 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_167 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_168 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_181 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_171 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_172 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_173 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_174 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_175 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_176 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_177 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_178 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
         const double xi_179 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_182 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_183 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_184 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_164 = 1.0*xi_159*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_165 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_166 = xi_159*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_161 + xi_164 + xi_165 + xi_166 + xi_170 + xi_180;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_167 + xi_168 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_181;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (coarse_level)) - 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
      {
         const double xi_71 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_87 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_79 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_77 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_88 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_83 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_84 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_85 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
         const double xi_86 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_89 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_90 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_91 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_74 = 1.0*xi_66*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_75 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_76 = xi_66*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_80 = 1.0*xi_69*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_81 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_82 = xi_69*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_71 + xi_74 + xi_75 + xi_76 + xi_79 + xi_87;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_77 + xi_79 + xi_80 + xi_81 + xi_82 + xi_88;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91;
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)) - 1; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_223 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_224 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_242 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_226 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_227 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_235 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_229 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_230 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_231 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_232 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_233 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_243 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_236 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_237 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_238 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_239 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_240 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_241 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_244 = 1.0*xi_221*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_245 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_246 = xi_221*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_223 + xi_224 + xi_226 + xi_227 + xi_229 + xi_230 + xi_231 + xi_235 + xi_242;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_232 + xi_233 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_243;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)) - 1; ctr_1 += 1)
      {
         const double xi_3 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_4 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_25 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_6 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_7 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_15 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_9 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_10 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_12 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_13 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_26 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_16 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_17 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_18 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_19 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_20 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_21 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_22 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_23 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
         const double xi_24 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_27 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_28 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_29 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_15 + xi_25 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (coarse_level)) - 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
      {
         const double xi_192 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_193 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_211 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_195 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_196 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_203 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_198 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_199 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_200 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_201 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_212 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_207 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_208 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_209 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
         const double xi_210 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_213 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_214 = 0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_215 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_204 = 1.0*xi_190*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_205 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_206 = xi_190*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_192 + xi_193 + xi_195 + xi_196 + xi_198 + xi_199 + xi_200 + xi_203 + xi_211;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_201 + xi_203 + xi_204 + xi_205 + xi_206 + xi_212;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_207 + xi_208 + xi_209 + xi_210 + xi_211 + xi_212 + xi_213 + xi_214 + xi_215;
      }
   }
   for (int ctr_2 = (1 << (coarse_level)) - 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_102 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_103 = 0.25*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_118 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_105 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_106 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_113 = 0.5*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_108 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_109 = 0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_110 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_111 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_119 = 0.5*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_117 = 0.25*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_114 = 1.0*xi_97*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_115 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_116 = xi_97*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_120 = 1.0*xi_100*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_121 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_122 = xi_100*0.75*_data_edgeFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_102 + xi_103 + xi_105 + xi_106 + xi_108 + xi_109 + xi_110 + xi_113 + xi_118;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_111 + xi_113 + xi_114 + xi_115 + xi_116 + xi_119;
         _data_edgeCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level)) + 1)*(1 << (coarse_level))) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122;
      }
   }
   {
      
   }
}


void restrict_2D_macroface_P2_update_edgedofs(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
    switch( coarse_level )
    {
    case 2:
        restrict_2D_macroface_P2_update_edgedofs_level_2(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 3:
        restrict_2D_macroface_P2_update_edgedofs_level_3(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 4:
        restrict_2D_macroface_P2_update_edgedofs_level_4(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 5:
        restrict_2D_macroface_P2_update_edgedofs_level_5(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 6:
        restrict_2D_macroface_P2_update_edgedofs_level_6(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 7:
        restrict_2D_macroface_P2_update_edgedofs_level_7(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 8:
        restrict_2D_macroface_P2_update_edgedofs_level_8(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 9:
        restrict_2D_macroface_P2_update_edgedofs_level_9(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 10:
        restrict_2D_macroface_P2_update_edgedofs_level_10(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 11:
        restrict_2D_macroface_P2_update_edgedofs_level_11(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 12:
        restrict_2D_macroface_P2_update_edgedofs_level_12(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 13:
        restrict_2D_macroface_P2_update_edgedofs_level_13(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    case 14:
        restrict_2D_macroface_P2_update_edgedofs_level_14(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    default:
        restrict_2D_macroface_P2_update_edgedofs_level_any(_data_edgeCoarseDst, _data_edgeFineSrc, _data_vertexFineSrc, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
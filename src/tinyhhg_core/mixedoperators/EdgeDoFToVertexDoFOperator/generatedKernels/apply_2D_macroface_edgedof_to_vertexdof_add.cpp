
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToVertexMacroFace2D.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_2(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 6] + xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_3(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 10] + xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_4(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 18] + xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_5(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 34] + xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_6(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 66] + xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_7(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 130] + xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_8(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 258] + xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_9(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 514] + xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_10(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026] + xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_11(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050] + xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_12(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098] + xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_13(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194] + xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_14(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386] + xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_add_level_any(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst, int64_t level)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[0];
   const double xi_2 = _data_edgeToVertexFaceStencil[10];
   const double xi_3 = _data_edgeToVertexFaceStencil[11];
   const double xi_4 = _data_edgeToVertexFaceStencil[9];
   const double xi_5 = _data_edgeToVertexFaceStencil[3];
   const double xi_6 = _data_edgeToVertexFaceStencil[4];
   const double xi_7 = _data_edgeToVertexFaceStencil[2];
   const double xi_8 = _data_edgeToVertexFaceStencil[6];
   const double xi_9 = _data_edgeToVertexFaceStencil[7];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < (1 << level); ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
      {
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_5*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + _data_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void apply_2D_macroface_edgedof_to_vertexdof_add(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_2(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 3:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_3(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 4:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_4(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 5:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_5(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 6:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_6(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 7:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_7(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 8:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_8(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 9:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_9(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 10:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_10(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 11:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_11(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 12:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_12(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 13:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_13(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 14:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_14(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    default:
        apply_2D_macroface_edgedof_to_vertexdof_add_level_any(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
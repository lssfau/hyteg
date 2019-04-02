
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToVertexMacroFace2D.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_2(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 6];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         _data_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_3(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 10];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         _data_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_4(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 18];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         _data_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_5(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 34];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         _data_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_6(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 66];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         _data_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_7(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 130];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         _data_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_8(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 258];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         _data_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_9(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 514];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         _data_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_10(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         _data_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_11(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         _data_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_12(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         _data_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_13(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         _data_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_14(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst)
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
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         _data_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst, int64_t level)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_15 = xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_19 = xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_20 = xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_21 = xi_5*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_22 = xi_6*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_23 = xi_7*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_24 = xi_8*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_16 = xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_17 = xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25;
      }
   }
}


void apply_2D_macroface_edgedof_to_vertexdof_replace(double * _data_edgeFaceSrc, double * const _data_edgeToVertexFaceStencil, double * _data_p1FaceDst, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_2(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 3:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_3(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 4:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_4(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 5:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_5(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 6:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_6(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 7:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_7(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 8:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_8(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 9:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_9(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 10:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_10(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 11:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_11(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 12:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_12(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 13:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_13(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    case 14:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_14(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst);
        break;
    default:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(_data_edgeFaceSrc, _data_edgeToVertexFaceStencil, _data_p1FaceDst, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
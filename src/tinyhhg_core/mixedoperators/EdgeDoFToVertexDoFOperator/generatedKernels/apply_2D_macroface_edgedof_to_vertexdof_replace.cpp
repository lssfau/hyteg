
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernels.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_2(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 4; ctr_1 < 5; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 6] + xi_10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + xi_2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 4; ctr_1 < -ctr_2 + 5; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 4; ctr_2 < 5; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_3(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 8; ctr_1 < 9; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 10] + xi_10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + xi_2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 8; ctr_1 < -ctr_2 + 9; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 8; ctr_2 < 9; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_4(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 16; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 16; ctr_1 < 17; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 18] + xi_10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + xi_2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 16; ctr_1 < -ctr_2 + 17; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 16; ctr_2 < 17; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_5(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 32; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 32; ctr_1 < 33; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 34] + xi_10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + xi_2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 32; ctr_1 < -ctr_2 + 33; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 32; ctr_2 < 33; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_6(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 64; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 64; ctr_1 < 65; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 66] + xi_10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + xi_2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 64; ctr_1 < -ctr_2 + 65; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 64; ctr_2 < 65; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_7(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 128; ctr_1 < 129; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 130] + xi_10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + xi_2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 128; ctr_1 < -ctr_2 + 129; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 128; ctr_2 < 129; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_8(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 256; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 256; ctr_1 < 257; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 258] + xi_10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + xi_2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 256; ctr_1 < -ctr_2 + 257; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 256; ctr_2 < 257; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_9(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 512; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 512; ctr_1 < 513; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 514] + xi_10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + xi_2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 512; ctr_1 < -ctr_2 + 513; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 512; ctr_2 < 513; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_10(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1024; ctr_1 < 1025; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026] + xi_10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + xi_2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 1024; ctr_1 < -ctr_2 + 1025; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1024; ctr_2 < 1025; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_11(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 2048; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 2048; ctr_1 < 2049; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050] + xi_10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + xi_2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 2048; ctr_1 < -ctr_2 + 2049; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 2048; ctr_2 < 2049; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_12(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 4096; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 4096; ctr_1 < 4097; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098] + xi_10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + xi_2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 4096; ctr_1 < -ctr_2 + 4097; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 4096; ctr_2 < 4097; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_13(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 8192; ctr_1 < 8193; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194] + xi_10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + xi_2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 8192; ctr_1 < -ctr_2 + 8193; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 8192; ctr_2 < 8193; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_14(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < 16384; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 16384; ctr_1 < 16385; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386] + xi_10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + xi_2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 16384; ctr_1 < -ctr_2 + 16385; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 16384; ctr_2 < 16385; ctr_2 += 1)
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
{
   const double xi_0 = fd_edgeToVertexFaceStencil[10];
   const double xi_1 = fd_edgeToVertexFaceStencil[1];
   const double xi_2 = fd_edgeToVertexFaceStencil[0];
   const double xi_3 = fd_edgeToVertexFaceStencil[3];
   const double xi_4 = fd_edgeToVertexFaceStencil[4];
   const double xi_5 = fd_edgeToVertexFaceStencil[11];
   const double xi_6 = fd_edgeToVertexFaceStencil[7];
   const double xi_7 = fd_edgeToVertexFaceStencil[2];
   const double xi_8 = fd_edgeToVertexFaceStencil[9];
   const double xi_9 = fd_edgeToVertexFaceStencil[8];
   const double xi_10 = fd_edgeToVertexFaceStencil[6];
   const double xi_11 = fd_edgeToVertexFaceStencil[5];
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = 1; ctr_1 < (1 << level); ctr_1 += 1)
   {
      
   }
   for (int ctr_1 = (1 << level); ctr_1 < (1 << level) + 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = 1; ctr_2 < (1 << level); ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
      {
         fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_11*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + xi_2*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_3*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + xi_4*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + (1 << level); ctr_1 < -ctr_2 + (1 << level) + 1; ctr_1 += 1)
      {
         
      }
   }
   for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
   {
      
   }
   for (int ctr_2 = (1 << level); ctr_2 < (1 << level) + 1; ctr_2 += 1)
   {
      
   }
}


void apply_2D_macroface_edgedof_to_vertexdof_replace(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_2(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 3:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_3(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 4:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_4(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 5:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_5(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 6:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_6(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 7:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_7(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 8:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_8(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 9:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_9(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 10:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_10(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 11:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_11(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 12:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_12(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 13:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_13(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    case 14:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_14(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
        break;
    default:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
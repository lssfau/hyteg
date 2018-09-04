
#include "generatedKernels.hpp"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_2(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (20 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (20 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(20 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 6] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 5];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(20 / 2) + 3] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 3] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 4] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 9] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 8];
    }
    for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[6*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6];
        fd_edgeFaceDst[5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6];
        fd_edgeFaceDst[ctr_1 + 5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 5];
      }
      {
        fd_edgeFaceDst[4*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 3] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 3] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[5*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[5*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9];
        fd_edgeFaceDst[4*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 3] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 3] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[5*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[5*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8];
      }
    }
    {
      fd_edgeFaceDst[-(12 / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 24];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (12 / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 24];
      }
      fd_edgeFaceDst[-(12 / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 24];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_3(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (72 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (72 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(72 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 9];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(72 / 2) + 7] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 7] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 17] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 16];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 9] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 11] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 10] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 20];
          fd_edgeFaceDst[-(2 / 2) + (72 / 2) + 9] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 21] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 20] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 10];
        }
        for (int ctr_1 = 1; ctr_1 < 6; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 9] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 11] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 20];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (72 / 2) + 9] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 21] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 20] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(72 / 2) + 9] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 11] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 20] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 19];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 17] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 16] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 26];
          fd_edgeFaceDst[-(2 / 2) + 2*(72 / 2) + 15] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 16] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 17] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 26] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 25];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 6; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10];
          fd_edgeFaceDst[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10];
          fd_edgeFaceDst[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9];
        }
        {
          fd_edgeFaceDst[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17];
          fd_edgeFaceDst[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16];
        }
      }
      {
        {
          fd_edgeFaceDst[-(42 / 2) + 54] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(42 / 2) + 60] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 51] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(56 / 2) + 70];
          fd_edgeFaceDst[-(42 / 2) + (72 / 2) + 54] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(56 / 2) + 71] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(56 / 2) + 70] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(42 / 2) + 60];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (42 / 2) + 54] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (42 / 2) + 60] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 51] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (56 / 2) + 70];
        }
        {
          fd_edgeFaceDst[-(42 / 2) + 55] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(42 / 2) + 62] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 52] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(56 / 2) + 71];
          fd_edgeFaceDst[-(42 / 2) + 2*(72 / 2) + 55] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(42 / 2) + 62] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(56 / 2) + 71] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(56 / 2) + 70];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (56 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (56 / 2) + 71] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (56 / 2) + 70] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (72 / 2) + 80];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_4(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (272 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (272 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(272 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 17];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(272 / 2) + 15] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 15] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 16] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 32];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 17] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 36];
          fd_edgeFaceDst[-(2 / 2) + (272 / 2) + 17] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 18];
        }
        for (int ctr_1 = 1; ctr_1 < 14; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 17] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (272 / 2) + 17] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(272 / 2) + 17] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 35];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 31] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 15] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 50];
          fd_edgeFaceDst[-(2 / 2) + 2*(272 / 2) + 31] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 49];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 34] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 54];
            fd_edgeFaceDst[-(6 / 2) + (272 / 2) + 34] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 54] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 36];
          }
          for (int ctr_1 = 1; ctr_1 < 13; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 34] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (272 / 2) + 34] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(272 / 2) + 34] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 53];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 47] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 67];
            fd_edgeFaceDst[-(6 / 2) + 2*(272 / 2) + 47] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 67] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 66];
          }
        }
        for (int ctr_2 = 3; ctr_2 < 13; ctr_2 += 1)
        {
          {
            fd_edgeFaceDst[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18];
            fd_edgeFaceDst[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          }
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18];
            fd_edgeFaceDst[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17];
          }
          {
            fd_edgeFaceDst[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33];
            fd_edgeFaceDst[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32];
          }
        }
        {
          {
            fd_edgeFaceDst[-(182 / 2) + 221] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 234] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 217] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 252];
            fd_edgeFaceDst[-(182 / 2) + (272 / 2) + 221] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 234];
          }
          for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (182 / 2) + 221] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 217] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 252];
            fd_edgeFaceDst[ctr_1 - (182 / 2) + (272 / 2) + 221] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 235] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 252] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234];
            fd_edgeFaceDst[ctr_1 - (182 / 2) + 2*(272 / 2) + 221] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 235] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 252] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 251];
          }
          {
            fd_edgeFaceDst[-(182 / 2) + 223] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 219] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 254];
            fd_edgeFaceDst[-(182 / 2) + 2*(272 / 2) + 223] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 253];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(210 / 2) + 238] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 270];
          fd_edgeFaceDst[-(210 / 2) + (272 / 2) + 238] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(240 / 2) + 271] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(240 / 2) + 270] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 252];
        }
        {
          fd_edgeFaceDst[-(210 / 2) + 239] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 271];
          fd_edgeFaceDst[-(210 / 2) + 2*(272 / 2) + 239] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(240 / 2) + 271] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 270];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (240 / 2) + 255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (240 / 2) + 271] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (240 / 2) + 270] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (272 / 2) + 288];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_5(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1056 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 33];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(1056 / 2) + 31] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 31] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 32] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 65] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 64];
    }
    for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34];
        fd_edgeFaceDst[33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34];
        fd_edgeFaceDst[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33];
      }
      {
        fd_edgeFaceDst[32*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 32] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65];
        fd_edgeFaceDst[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 32] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 64];
      }
    }
    {
      fd_edgeFaceDst[-(992 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(992 / 2) + 1055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(992 / 2) + 1054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1056 / 2) + 1088];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (992 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (992 / 2) + 1055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (992 / 2) + 1054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1056 / 2) + 1088];
      }
      fd_edgeFaceDst[-(992 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(992 / 2) + 1055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(992 / 2) + 1054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1056 / 2) + 1088];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_6(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4160 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4160 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4160 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 65];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(4160 / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 63] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 64] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 128];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 65] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 132];
          fd_edgeFaceDst[-(2 / 2) + (4160 / 2) + 65] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 132] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 66];
        }
        for (int ctr_1 = 1; ctr_1 < 62; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 65] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (4160 / 2) + 65] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(4160 / 2) + 65] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 131];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 63] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 194];
          fd_edgeFaceDst[-(2 / 2) + 2*(4160 / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 193];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 62; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66];
          fd_edgeFaceDst[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66];
          fd_edgeFaceDst[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65];
        }
        {
          fd_edgeFaceDst[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 64] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129];
          fd_edgeFaceDst[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 64] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128];
        }
      }
      {
        {
          fd_edgeFaceDst[-(3906 / 2) + 4030] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4092] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3782 / 2) + 4027] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4158];
          fd_edgeFaceDst[-(3906 / 2) + (4160 / 2) + 4030] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4032 / 2) + 4159] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4032 / 2) + 4158] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(3906 / 2) + 4092];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (3906 / 2) + 4030] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4092] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4027] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4032 / 2) + 4158];
        }
        {
          fd_edgeFaceDst[-(3906 / 2) + 4031] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3782 / 2) + 4028] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4159];
          fd_edgeFaceDst[-(3906 / 2) + 2*(4160 / 2) + 4031] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4032 / 2) + 4159] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4158];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (4032 / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4032 / 2) + 4159] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4032 / 2) + 4158] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4160 / 2) + 4224];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_7(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16512 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16512 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16512 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 129];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(16512 / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 127] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 257] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 256];
    }
    for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130];
        fd_edgeFaceDst[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130];
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129];
      }
      {
        fd_edgeFaceDst[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 128] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257];
        fd_edgeFaceDst[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256];
      }
    }
    {
      fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16256 / 2) + 16511] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16256 / 2) + 16510] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16512 / 2) + 16640];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (16256 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16256 / 2) + 16511] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16256 / 2) + 16510] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16512 / 2) + 16640];
      }
      fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16256 / 2) + 16511] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16256 / 2) + 16510] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16512 / 2) + 16640];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_8(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (65792 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (65792 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(65792 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 257];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(65792 / 2) + 255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 255] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 256] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 513] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 512];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 257] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 259] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 258] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 516];
          fd_edgeFaceDst[-(2 / 2) + (65792 / 2) + 257] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 517] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 516] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 258];
        }
        for (int ctr_1 = 1; ctr_1 < 254; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 257] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (65792 / 2) + 257] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 517] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(65792 / 2) + 257] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 515];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 511] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 513] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 512] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 255] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 770];
          fd_edgeFaceDst[-(2 / 2) + 2*(65792 / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 512] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 513] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 770] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 769];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 514] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 517] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 516] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 259] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 774];
            fd_edgeFaceDst[-(6 / 2) + (65792 / 2) + 514] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 517] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 775] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 774] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 516];
          }
          for (int ctr_1 = 1; ctr_1 < 253; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 514] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 517] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (65792 / 2) + 514] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 517] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 775] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(65792 / 2) + 514] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 516] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 517] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 773];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 767] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 770] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 769] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 512] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 1027];
            fd_edgeFaceDst[-(6 / 2) + 2*(65792 / 2) + 767] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 769] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 770] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 1027] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 1026];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 771] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 775] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 774] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 517] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 1032];
              fd_edgeFaceDst[-(12 / 2) + (65792 / 2) + 771] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 775] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 1033] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 1032] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 774];
            }
            for (int ctr_1 = 1; ctr_1 < 252; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 771] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 775] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 517] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 1032];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (65792 / 2) + 771] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 775] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 1033] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 1032] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(65792 / 2) + 771] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 774] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 775] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 1032] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 1031];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 1027] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 1026] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 769] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 1284];
              fd_edgeFaceDst[-(12 / 2) + 2*(65792 / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 1027] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 1284] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 1283];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 252; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258];
              fd_edgeFaceDst[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258];
              fd_edgeFaceDst[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257];
            }
            {
              fd_edgeFaceDst[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 256] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513];
              fd_edgeFaceDst[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 256] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512];
            }
          }
          {
            {
              fd_edgeFaceDst[-(63756 / 2) + 64764] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(63756 / 2) + 65017] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(63756 / 2) + 65016] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(63252 / 2) + 64759] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(64262 / 2) + 65274];
              fd_edgeFaceDst[-(63756 / 2) + (65792 / 2) + 64764] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(63756 / 2) + 65017] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(64262 / 2) + 65274] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(63756 / 2) + 65016];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (63756 / 2) + 64764] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65017] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65016] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (63252 / 2) + 64759] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (64262 / 2) + 65274];
              fd_edgeFaceDst[ctr_1 - (63756 / 2) + (65792 / 2) + 64764] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65017] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (64262 / 2) + 65275] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (64262 / 2) + 65274] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65016];
              fd_edgeFaceDst[ctr_1 - (63756 / 2) + 2*(65792 / 2) + 64764] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65016] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (63756 / 2) + 65017] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (64262 / 2) + 65274] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (64262 / 2) + 65273];
            }
            {
              fd_edgeFaceDst[-(63756 / 2) + 64767] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(63756 / 2) + 65020] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(63756 / 2) + 65019] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(63252 / 2) + 64762] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(64262 / 2) + 65277];
              fd_edgeFaceDst[-(63756 / 2) + 2*(65792 / 2) + 64767] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(63756 / 2) + 65019] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(63756 / 2) + 65020] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(64262 / 2) + 65277] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(64262 / 2) + 65276];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(64262 / 2) + 65021] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65274] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(63756 / 2) + 65017] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65532];
            fd_edgeFaceDst[-(64262 / 2) + (65792 / 2) + 65021] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65532] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(64262 / 2) + 65274];
          }
          {
            fd_edgeFaceDst[-(64262 / 2) + 65022] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(63756 / 2) + 65018] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65533];
            fd_edgeFaceDst[-(64262 / 2) + (65792 / 2) + 65022] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(64770 / 2) + 65534] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(64262 / 2) + 65275];
            fd_edgeFaceDst[-(64262 / 2) + 2*(65792 / 2) + 65022] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65532];
          }
          {
            fd_edgeFaceDst[-(64262 / 2) + 65023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65277] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(63756 / 2) + 65019] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65534];
            fd_edgeFaceDst[-(64262 / 2) + 2*(65792 / 2) + 65023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(64262 / 2) + 65277] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65534] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65533];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(64770 / 2) + 65278] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(64770 / 2) + 65532] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(64262 / 2) + 65275] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(65280 / 2) + 65790];
          fd_edgeFaceDst[-(64770 / 2) + (65792 / 2) + 65278] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(65280 / 2) + 65791] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(65280 / 2) + 65790] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(64770 / 2) + 65532];
        }
        {
          fd_edgeFaceDst[-(64770 / 2) + 65279] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(64770 / 2) + 65534] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(64262 / 2) + 65276] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(65280 / 2) + 65791];
          fd_edgeFaceDst[-(64770 / 2) + 2*(65792 / 2) + 65279] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(64770 / 2) + 65534] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(65280 / 2) + 65791] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(65280 / 2) + 65790];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (65280 / 2) + 65535] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (65280 / 2) + 65791] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (65280 / 2) + 65790] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (65792 / 2) + 66048];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_9(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (262656 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (262656 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(262656 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 513];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(262656 / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 511] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 512] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 1025] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 1024];
    }
    for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514];
        fd_edgeFaceDst[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514];
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513];
      }
      {
        fd_edgeFaceDst[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 512] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025];
        fd_edgeFaceDst[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 512] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024];
      }
    }
    {
      fd_edgeFaceDst[-(261632 / 2) + 262143] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(261632 / 2) + 262655] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(261632 / 2) + 262654] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(262656 / 2) + 263168];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (261632 / 2) + 262143] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (261632 / 2) + 262655] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (261632 / 2) + 262654] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (262656 / 2) + 263168];
      }
      fd_edgeFaceDst[-(261632 / 2) + 262143] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(261632 / 2) + 262655] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(261632 / 2) + 262654] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(262656 / 2) + 263168];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_10(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1049600 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1049600 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1049600 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1025];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(1049600 / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1023] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 1024] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2049] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 2048];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 1025] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 1027] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 1026] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 2052];
          fd_edgeFaceDst[-(2 / 2) + (1049600 / 2) + 1025] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 2053] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 2052] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 1026];
        }
        for (int ctr_1 = 1; ctr_1 < 1022; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 1025] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2053] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2051];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 2047] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 2049] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 2048] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1023] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 3074];
          fd_edgeFaceDst[-(2 / 2) + 2*(1049600 / 2) + 2047] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 2048] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 2049] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 3074] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 3073];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 2050] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 2053] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 2052] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 1027] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 3078];
            fd_edgeFaceDst[-(6 / 2) + (1049600 / 2) + 2050] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 2053] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 3079] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 3078] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 2052];
          }
          for (int ctr_1 = 1; ctr_1 < 1021; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2050] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2053] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2053] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3079] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2052] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2053] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3077];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 3071] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 3074] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 3073] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2048] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 4099];
            fd_edgeFaceDst[-(6 / 2) + 2*(1049600 / 2) + 3071] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 3073] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 3074] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 4099] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 4098];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 3075] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 3079] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 3078] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 2053] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 4104];
              fd_edgeFaceDst[-(12 / 2) + (1049600 / 2) + 3075] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 3079] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 4105] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 4104] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 3078];
            }
            for (int ctr_1 = 1; ctr_1 < 1020; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 3075] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3079] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 2053] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 4104];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (1049600 / 2) + 3075] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3079] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 4105] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 4104] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3075] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3078] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 3079] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 4104] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 4103];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 4099] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 4098] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 3073] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 5124];
              fd_edgeFaceDst[-(12 / 2) + 2*(1049600 / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 4099] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 5124] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 5123];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 1020; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026];
              fd_edgeFaceDst[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026];
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025];
            }
            {
              fd_edgeFaceDst[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1024] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049];
              fd_edgeFaceDst[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1024] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2048];
            }
          }
          {
            {
              fd_edgeFaceDst[-(1041420 / 2) + 1045500] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1041420 / 2) + 1046521] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1041420 / 2) + 1046520] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1039380 / 2) + 1045495] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1043462 / 2) + 1047546];
              fd_edgeFaceDst[-(1041420 / 2) + (1049600 / 2) + 1045500] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(1041420 / 2) + 1046521] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(1043462 / 2) + 1047546] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(1041420 / 2) + 1046520];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + 1045500] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046521] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046520] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (1039380 / 2) + 1045495] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1043462 / 2) + 1047546];
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046521] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (1043462 / 2) + 1047547] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (1043462 / 2) + 1047546] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046520];
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045500] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046520] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (1041420 / 2) + 1046521] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (1043462 / 2) + 1047546] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1043462 / 2) + 1047545];
            }
            {
              fd_edgeFaceDst[-(1041420 / 2) + 1045503] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1041420 / 2) + 1046524] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1041420 / 2) + 1046523] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1039380 / 2) + 1045498] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1043462 / 2) + 1047549];
              fd_edgeFaceDst[-(1041420 / 2) + 2*(1049600 / 2) + 1045503] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(1041420 / 2) + 1046523] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(1041420 / 2) + 1046524] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(1043462 / 2) + 1047549] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(1043462 / 2) + 1047548];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046525] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047546] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1041420 / 2) + 1046521] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048572];
            fd_edgeFaceDst[-(1043462 / 2) + (1049600 / 2) + 1046525] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048572] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(1043462 / 2) + 1047546];
          }
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046526] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1041420 / 2) + 1046522] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048573];
            fd_edgeFaceDst[-(1043462 / 2) + (1049600 / 2) + 1046526] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(1045506 / 2) + 1048574] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(1043462 / 2) + 1047547];
            fd_edgeFaceDst[-(1043462 / 2) + 2*(1049600 / 2) + 1046526] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048572];
          }
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046527] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047549] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1041420 / 2) + 1046523] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048574];
            fd_edgeFaceDst[-(1043462 / 2) + 2*(1049600 / 2) + 1046527] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(1043462 / 2) + 1047549] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048574] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048573];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(1045506 / 2) + 1047550] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1045506 / 2) + 1048572] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1043462 / 2) + 1047547] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1047552 / 2) + 1049598];
          fd_edgeFaceDst[-(1045506 / 2) + (1049600 / 2) + 1047550] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(1047552 / 2) + 1049599] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(1047552 / 2) + 1049598] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(1045506 / 2) + 1048572];
        }
        {
          fd_edgeFaceDst[-(1045506 / 2) + 1047551] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1045506 / 2) + 1048574] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1043462 / 2) + 1047548] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1047552 / 2) + 1049599];
          fd_edgeFaceDst[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(1045506 / 2) + 1048574] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(1047552 / 2) + 1049599] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(1047552 / 2) + 1049598];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (1047552 / 2) + 1048575] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (1047552 / 2) + 1049599] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (1047552 / 2) + 1049598] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1049600 / 2) + 1050624];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_11(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4196352 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4196352 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4196352 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2049];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(4196352 / 2) + 2047] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 2047] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 2048] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4097] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 4096];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 2049] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 2051] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 2050] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 4100];
          fd_edgeFaceDst[-(2 / 2) + (4196352 / 2) + 2049] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 4101] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 4100] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 2050];
        }
        for (int ctr_1 = 1; ctr_1 < 2046; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2049] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (4196352 / 2) + 2049] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4101] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(4196352 / 2) + 2049] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4099];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 4097] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 4096] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 2047] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 6146];
          fd_edgeFaceDst[-(2 / 2) + 2*(4196352 / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 4096] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 4097] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 6146] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 6145];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 4098] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 4101] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 4100] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2051] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 6150];
            fd_edgeFaceDst[-(6 / 2) + (4196352 / 2) + 4098] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 4101] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 6151] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 6150] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 4100];
          }
          for (int ctr_1 = 1; ctr_1 < 2045; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 4098] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4101] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (4196352 / 2) + 4098] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4101] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6151] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(4196352 / 2) + 4098] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4100] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4101] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6149];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 6143] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 6146] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 6145] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4096] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 8195];
            fd_edgeFaceDst[-(6 / 2) + 2*(4196352 / 2) + 6143] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 6145] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 6146] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 8195] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 8194];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 6147] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 6151] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 6150] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 4101] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 8200];
              fd_edgeFaceDst[-(12 / 2) + (4196352 / 2) + 6147] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 6151] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 8201] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 8200] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 6150];
            }
            for (int ctr_1 = 1; ctr_1 < 2044; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 6147] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6151] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 4101] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 8200];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (4196352 / 2) + 6147] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6151] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 8201] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 8200] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(4196352 / 2) + 6147] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6150] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 6151] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 8200] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 8199];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 8191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 8195] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 8194] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 6145] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 10244];
              fd_edgeFaceDst[-(12 / 2) + 2*(4196352 / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 8195] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 10244] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 10243];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 2044; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050];
              fd_edgeFaceDst[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050];
              fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049];
            }
            {
              fd_edgeFaceDst[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2048] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097];
              fd_edgeFaceDst[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2048] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096];
            }
          }
          {
            {
              fd_edgeFaceDst[-(4179980 / 2) + 4188156] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4179980 / 2) + 4190201] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4179980 / 2) + 4190200] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4175892 / 2) + 4188151] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4184070 / 2) + 4192250];
              fd_edgeFaceDst[-(4179980 / 2) + (4196352 / 2) + 4188156] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(4179980 / 2) + 4190201] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4184070 / 2) + 4192250] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(4179980 / 2) + 4190200];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (4179980 / 2) + 4188156] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190201] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190200] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (4175892 / 2) + 4188151] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4184070 / 2) + 4192250];
              fd_edgeFaceDst[ctr_1 - (4179980 / 2) + (4196352 / 2) + 4188156] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190201] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4184070 / 2) + 4192251] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (4184070 / 2) + 4192250] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190200];
              fd_edgeFaceDst[ctr_1 - (4179980 / 2) + 2*(4196352 / 2) + 4188156] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190200] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4179980 / 2) + 4190201] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (4184070 / 2) + 4192250] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4184070 / 2) + 4192249];
            }
            {
              fd_edgeFaceDst[-(4179980 / 2) + 4188159] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4179980 / 2) + 4190204] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4179980 / 2) + 4190203] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4175892 / 2) + 4188154] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4184070 / 2) + 4192253];
              fd_edgeFaceDst[-(4179980 / 2) + 2*(4196352 / 2) + 4188159] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(4179980 / 2) + 4190203] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(4179980 / 2) + 4190204] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4184070 / 2) + 4192253] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4184070 / 2) + 4192252];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(4184070 / 2) + 4190205] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192250] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4179980 / 2) + 4190201] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194300];
            fd_edgeFaceDst[-(4184070 / 2) + (4196352 / 2) + 4190205] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194300] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(4184070 / 2) + 4192250];
          }
          {
            fd_edgeFaceDst[-(4184070 / 2) + 4190206] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4179980 / 2) + 4190202] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194301];
            fd_edgeFaceDst[-(4184070 / 2) + (4196352 / 2) + 4190206] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4188162 / 2) + 4194302] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(4184070 / 2) + 4192251];
            fd_edgeFaceDst[-(4184070 / 2) + 2*(4196352 / 2) + 4190206] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194300];
          }
          {
            fd_edgeFaceDst[-(4184070 / 2) + 4190207] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192253] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4179980 / 2) + 4190203] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194302];
            fd_edgeFaceDst[-(4184070 / 2) + 2*(4196352 / 2) + 4190207] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(4184070 / 2) + 4192253] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194302] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194301];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(4188162 / 2) + 4192254] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4188162 / 2) + 4194300] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4184070 / 2) + 4192251] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4192256 / 2) + 4196350];
          fd_edgeFaceDst[-(4188162 / 2) + (4196352 / 2) + 4192254] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4192256 / 2) + 4196351] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4192256 / 2) + 4196350] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(4188162 / 2) + 4194300];
        }
        {
          fd_edgeFaceDst[-(4188162 / 2) + 4192255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4188162 / 2) + 4194302] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4184070 / 2) + 4192252] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4192256 / 2) + 4196351];
          fd_edgeFaceDst[-(4188162 / 2) + 2*(4196352 / 2) + 4192255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(4188162 / 2) + 4194302] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4192256 / 2) + 4196351] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4192256 / 2) + 4196350];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (4192256 / 2) + 4194303] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4192256 / 2) + 4196351] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4192256 / 2) + 4196350] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4196352 / 2) + 4198400];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_12(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16781312 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16781312 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16781312 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4097];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(16781312 / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 4095] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 4096] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8193] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 8192];
    }
    for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098];
        fd_edgeFaceDst[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098];
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097];
      }
      {
        fd_edgeFaceDst[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4096] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193];
        fd_edgeFaceDst[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4096] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192];
      }
    }
    {
      fd_edgeFaceDst[-(16773120 / 2) + 16777215] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16773120 / 2) + 16781311] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16773120 / 2) + 16781310] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16781312 / 2) + 16785408];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (16773120 / 2) + 16777215] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16773120 / 2) + 16781311] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16773120 / 2) + 16781310] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16781312 / 2) + 16785408];
      }
      fd_edgeFaceDst[-(16773120 / 2) + 16777215] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16773120 / 2) + 16781311] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16773120 / 2) + 16781310] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16781312 / 2) + 16785408];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_13(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (67117056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (67117056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(67117056 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8193];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(67117056 / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 8191] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 8192] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 16385] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 16384];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 8193] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8195] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8194] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 16388];
          fd_edgeFaceDst[-(2 / 2) + (67117056 / 2) + 8193] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 16389] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 16388] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 8194];
        }
        for (int ctr_1 = 1; ctr_1 < 8190; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 8193] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16389] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16387];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 16385] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 16384] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 8191] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 24578];
          fd_edgeFaceDst[-(2 / 2) + 2*(67117056 / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 16384] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 16385] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 24578] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 24577];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 16386] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 16389] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 16388] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8195] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 24582];
            fd_edgeFaceDst[-(6 / 2) + (67117056 / 2) + 16386] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 16389] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 24583] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 24582] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 16388];
          }
          for (int ctr_1 = 1; ctr_1 < 8189; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 16386] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16389] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16389] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24583] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16388] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16389] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24581];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 24575] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 24578] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 24577] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 16384] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 32771];
            fd_edgeFaceDst[-(6 / 2) + 2*(67117056 / 2) + 24575] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 24577] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 24578] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 32771] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 32770];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 24579] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 24583] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 24582] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 16389] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 32776];
              fd_edgeFaceDst[-(12 / 2) + (67117056 / 2) + 24579] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 24583] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 32777] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 32776] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 24582];
            }
            for (int ctr_1 = 1; ctr_1 < 8188; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 24579] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24583] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 16389] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 32776];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (67117056 / 2) + 24579] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24583] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 32777] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 32776] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(67117056 / 2) + 24579] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24582] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 24583] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 32776] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 32775];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 32767] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 32771] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 32770] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 24577] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 40964];
              fd_edgeFaceDst[-(12 / 2) + 2*(67117056 / 2) + 32767] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 32770] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 32771] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 40964] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 40963];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 8188; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194];
              fd_edgeFaceDst[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194];
              fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193];
            }
            {
              fd_edgeFaceDst[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8192] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385];
              fd_edgeFaceDst[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8192] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384];
            }
          }
          {
            {
              fd_edgeFaceDst[-(67051532 / 2) + 67084284] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67051532 / 2) + 67092473] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67051532 / 2) + 67092472] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67035156 / 2) + 67084279] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67067910 / 2) + 67100666];
              fd_edgeFaceDst[-(67051532 / 2) + (67117056 / 2) + 67084284] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(67051532 / 2) + 67092473] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(67067910 / 2) + 67100666] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(67051532 / 2) + 67092472];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (67051532 / 2) + 67084284] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092473] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092472] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (67035156 / 2) + 67084279] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (67067910 / 2) + 67100666];
              fd_edgeFaceDst[ctr_1 - (67051532 / 2) + (67117056 / 2) + 67084284] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092473] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (67067910 / 2) + 67100667] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (67067910 / 2) + 67100666] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092472];
              fd_edgeFaceDst[ctr_1 - (67051532 / 2) + 2*(67117056 / 2) + 67084284] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092472] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (67051532 / 2) + 67092473] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (67067910 / 2) + 67100666] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (67067910 / 2) + 67100665];
            }
            {
              fd_edgeFaceDst[-(67051532 / 2) + 67084287] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67051532 / 2) + 67092476] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67051532 / 2) + 67092475] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67035156 / 2) + 67084282] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67067910 / 2) + 67100669];
              fd_edgeFaceDst[-(67051532 / 2) + 2*(67117056 / 2) + 67084287] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(67051532 / 2) + 67092475] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(67051532 / 2) + 67092476] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(67067910 / 2) + 67100669] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(67067910 / 2) + 67100668];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(67067910 / 2) + 67092477] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100666] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67051532 / 2) + 67092473] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108860];
            fd_edgeFaceDst[-(67067910 / 2) + (67117056 / 2) + 67092477] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108860] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(67067910 / 2) + 67100666];
          }
          {
            fd_edgeFaceDst[-(67067910 / 2) + 67092478] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67051532 / 2) + 67092474] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108861];
            fd_edgeFaceDst[-(67067910 / 2) + (67117056 / 2) + 67092478] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(67084290 / 2) + 67108862] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(67067910 / 2) + 67100667];
            fd_edgeFaceDst[-(67067910 / 2) + 2*(67117056 / 2) + 67092478] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108860];
          }
          {
            fd_edgeFaceDst[-(67067910 / 2) + 67092479] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100669] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67051532 / 2) + 67092475] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108862];
            fd_edgeFaceDst[-(67067910 / 2) + 2*(67117056 / 2) + 67092479] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(67067910 / 2) + 67100669] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108862] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108861];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(67084290 / 2) + 67100670] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67084290 / 2) + 67108860] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67067910 / 2) + 67100667] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67100672 / 2) + 67117054];
          fd_edgeFaceDst[-(67084290 / 2) + (67117056 / 2) + 67100670] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(67100672 / 2) + 67117055] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(67100672 / 2) + 67117054] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(67084290 / 2) + 67108860];
        }
        {
          fd_edgeFaceDst[-(67084290 / 2) + 67100671] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67084290 / 2) + 67108862] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67067910 / 2) + 67100668] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67100672 / 2) + 67117055];
          fd_edgeFaceDst[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(67084290 / 2) + 67108862] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(67100672 / 2) + 67117055] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(67100672 / 2) + 67117054];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (67100672 / 2) + 67108863] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (67100672 / 2) + 67117055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (67100672 / 2) + 67117054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (67117056 / 2) + 67125248];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_14(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (268451840 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)];
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (268451840 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(268451840 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16386] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16385];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(268451840 / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 16383] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 16384] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 32769] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 32768];
    }
    for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386];
        fd_edgeFaceDst[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386];
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385];
      }
      {
        fd_edgeFaceDst[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16384] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32769];
        fd_edgeFaceDst[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16384] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32769] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32768];
      }
    }
    {
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(268451840 / 2) + 268468224];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (268451840 / 2) + 268468224];
      }
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(268451840 / 2) + 268468224];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < (1 << level); ctr_2 += 1)
    for (int ctr_1 = 0; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      if (ctr_2 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)];
      }
      if (ctr_1 + ctr_2 < (1 << level) - 1)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      if (ctr_1 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1];
      }
    }
}




static void apply_2D_macroface_vertexdof_to_edgedof_replace(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  switch( level )
  {
#if 0
    case 2:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_2(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 3:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_3(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 4:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_4(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 5:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_5(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 6:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_6(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 7:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_7(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 8:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_8(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 9:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_9(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 10:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_10(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 11:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_11(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 12:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_12(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 13:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_13(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 14:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_14(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
#endif
    default:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil, level);
      break;
  }
}




static void apply_2D_macroface_vertexdof_to_edgedof_add_level_2(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (20 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (20 / 2)];
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (20 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 6] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (20 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(20 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 6] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(20 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(20 / 2) + 3] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 3] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 4] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 9] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_edgeFaceDst[-(0 / 2) + 2*(20 / 2) + 3];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 5] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 6] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 12] + fd_edgeFaceDst[-(2 / 2) + 5];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 5] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 12] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 6] + fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 5];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 6] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_edgeFaceDst[-(2 / 2) + 6];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 6] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 14] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 6];
          fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 6] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 12] + fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 6];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 7] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 9] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 3] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 14] + fd_edgeFaceDst[-(2 / 2) + 7];
          fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 7] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 9] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 14] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 7];
        }
      }
      {
        {
          fd_edgeFaceDst[-(6 / 2) + 10] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 12] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 18] + fd_edgeFaceDst[-(6 / 2) + 10];
          fd_edgeFaceDst[-(6 / 2) + (20 / 2) + 10] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 12] + fd_edgeFaceDst[-(6 / 2) + (20 / 2) + 10];
        }
        {
          fd_edgeFaceDst[-(6 / 2) + 11] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 14] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 19] + fd_edgeFaceDst[-(6 / 2) + 11];
          fd_edgeFaceDst[-(6 / 2) + 2*(20 / 2) + 11] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 13] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 14] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 19] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 18] + fd_edgeFaceDst[-(6 / 2) + 2*(20 / 2) + 11];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (12 / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 13] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 24] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 15];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_3(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (72 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (72 / 2)];
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (72 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (72 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(72 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 10] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 9] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(72 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(72 / 2) + 7] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 7] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 17] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 16] + fd_edgeFaceDst[-(0 / 2) + 2*(72 / 2) + 7];
    }
    for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_edgeFaceDst[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 11] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_edgeFaceDst[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7];
        fd_edgeFaceDst[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16] + fd_edgeFaceDst[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7];
      }
    }
    {
      fd_edgeFaceDst[-(56 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(56 / 2) + 71] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(56 / 2) + 70] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(72 / 2) + 80] + fd_edgeFaceDst[-(56 / 2) + 63];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (56 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (56 / 2) + 71] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (56 / 2) + 70] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (72 / 2) + 80] + fd_edgeFaceDst[ctr_1 - (56 / 2) + 63];
      }
      fd_edgeFaceDst[-(56 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(56 / 2) + 71] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(56 / 2) + 70] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(42 / 2) + 61] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(72 / 2) + 80] + fd_edgeFaceDst[-(56 / 2) + 63];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_4(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (272 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (272 / 2)];
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (272 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (272 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(272 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(272 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(272 / 2) + 15] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 15] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 16] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_edgeFaceDst[-(0 / 2) + 2*(272 / 2) + 15];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 17] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_edgeFaceDst[-(2 / 2) + 17];
          fd_edgeFaceDst[-(2 / 2) + (272 / 2) + 17] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 18] + fd_edgeFaceDst[-(2 / 2) + (272 / 2) + 17];
        }
        for (int ctr_1 = 1; ctr_1 < 14; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 17] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 17];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (272 / 2) + 17] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (272 / 2) + 17];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(272 / 2) + 17] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 35] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(272 / 2) + 17];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 31] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 15] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_edgeFaceDst[-(2 / 2) + 31];
          fd_edgeFaceDst[-(2 / 2) + 2*(272 / 2) + 31] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_edgeFaceDst[-(2 / 2) + 2*(272 / 2) + 31];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 34] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 54] + fd_edgeFaceDst[-(6 / 2) + 34];
            fd_edgeFaceDst[-(6 / 2) + (272 / 2) + 34] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 54] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 36] + fd_edgeFaceDst[-(6 / 2) + (272 / 2) + 34];
          }
          for (int ctr_1 = 1; ctr_1 < 13; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 34] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 19] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 34];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (272 / 2) + 34] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_edgeFaceDst[ctr_1 - (6 / 2) + (272 / 2) + 34];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(272 / 2) + 34] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 36] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 53] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(272 / 2) + 34];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 47] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 32] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 67] + fd_edgeFaceDst[-(6 / 2) + 47];
            fd_edgeFaceDst[-(6 / 2) + 2*(272 / 2) + 47] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 50] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 67] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 66] + fd_edgeFaceDst[-(6 / 2) + 2*(272 / 2) + 47];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 51] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 55] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 54] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 72] + fd_edgeFaceDst[-(12 / 2) + 51];
              fd_edgeFaceDst[-(12 / 2) + (272 / 2) + 51] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 73] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 72] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 54] + fd_edgeFaceDst[-(12 / 2) + (272 / 2) + 51];
            }
            for (int ctr_1 = 1; ctr_1 < 12; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 51] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 37] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 51];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (272 / 2) + 51] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 73] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_edgeFaceDst[ctr_1 - (12 / 2) + (272 / 2) + 51];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(272 / 2) + 51] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 54] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 71] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(272 / 2) + 51];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 49] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 84] + fd_edgeFaceDst[-(12 / 2) + 63];
              fd_edgeFaceDst[-(12 / 2) + 2*(272 / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 67] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 84] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 83] + fd_edgeFaceDst[-(12 / 2) + 2*(272 / 2) + 63];
            }
          }
          {
            {
              {
                fd_edgeFaceDst[-(20 / 2) + 68] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 73] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 72] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 55] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 90] + fd_edgeFaceDst[-(20 / 2) + 68];
                fd_edgeFaceDst[-(20 / 2) + (272 / 2) + 68] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 73] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(30 / 2) + 91] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 90] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 72] + fd_edgeFaceDst[-(20 / 2) + (272 / 2) + 68];
              }
              for (int ctr_1 = 1; ctr_1 < 11; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 68] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 73] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 55] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 90] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 68];
                fd_edgeFaceDst[ctr_1 - (20 / 2) + (272 / 2) + 68] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 73] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 91] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 90] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_edgeFaceDst[ctr_1 - (20 / 2) + (272 / 2) + 68];
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(272 / 2) + 68] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 72] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 73] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 90] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 89] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(272 / 2) + 68];
              }
              {
                fd_edgeFaceDst[-(20 / 2) + 79] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 84] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 83] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 101] + fd_edgeFaceDst[-(20 / 2) + 79];
                fd_edgeFaceDst[-(20 / 2) + 2*(272 / 2) + 79] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 83] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 84] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 101] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 100] + fd_edgeFaceDst[-(20 / 2) + 2*(272 / 2) + 79];
              }
            }
            for (int ctr_2 = 5; ctr_2 < 11; ctr_2 += 1)
            {
              {
                fd_edgeFaceDst[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_edgeFaceDst[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              }
              for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 19] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              }
              {
                fd_edgeFaceDst[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_edgeFaceDst[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15];
                fd_edgeFaceDst[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32] + fd_edgeFaceDst[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15];
              }
            }
            {
              {
                fd_edgeFaceDst[-(132 / 2) + 187] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(132 / 2) + 199] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(132 / 2) + 198] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(110 / 2) + 181] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(156 / 2) + 216] + fd_edgeFaceDst[-(132 / 2) + 187];
                fd_edgeFaceDst[-(132 / 2) + (272 / 2) + 187] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(132 / 2) + 199] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(156 / 2) + 217] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 216] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(132 / 2) + 198] + fd_edgeFaceDst[-(132 / 2) + (272 / 2) + 187];
              }
              for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 - (132 / 2) + 187] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 199] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 198] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (110 / 2) + 181] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_edgeFaceDst[ctr_1 - (132 / 2) + 187];
                fd_edgeFaceDst[ctr_1 - (132 / 2) + (272 / 2) + 187] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 199] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 217] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 198] + fd_edgeFaceDst[ctr_1 - (132 / 2) + (272 / 2) + 187];
                fd_edgeFaceDst[ctr_1 - (132 / 2) + 2*(272 / 2) + 187] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 198] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 199] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 215] + fd_edgeFaceDst[ctr_1 - (132 / 2) + 2*(272 / 2) + 187];
              }
              {
                fd_edgeFaceDst[-(132 / 2) + 191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(132 / 2) + 203] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(132 / 2) + 202] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(110 / 2) + 185] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(156 / 2) + 220] + fd_edgeFaceDst[-(132 / 2) + 191];
                fd_edgeFaceDst[-(132 / 2) + 2*(272 / 2) + 191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(132 / 2) + 202] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(132 / 2) + 203] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 220] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(156 / 2) + 219] + fd_edgeFaceDst[-(132 / 2) + 2*(272 / 2) + 191];
              }
            }
          }
          {
            {
              fd_edgeFaceDst[-(156 / 2) + 204] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(156 / 2) + 217] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(156 / 2) + 216] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(132 / 2) + 199] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 234] + fd_edgeFaceDst[-(156 / 2) + 204];
              fd_edgeFaceDst[-(156 / 2) + (272 / 2) + 204] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(156 / 2) + 217] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 234] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(156 / 2) + 216] + fd_edgeFaceDst[-(156 / 2) + (272 / 2) + 204];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (156 / 2) + 204] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 217] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (132 / 2) + 199] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234] + fd_edgeFaceDst[ctr_1 - (156 / 2) + 204];
              fd_edgeFaceDst[ctr_1 - (156 / 2) + (272 / 2) + 204] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 217] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 235] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_edgeFaceDst[ctr_1 - (156 / 2) + (272 / 2) + 204];
              fd_edgeFaceDst[ctr_1 - (156 / 2) + 2*(272 / 2) + 204] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 216] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (156 / 2) + 217] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 234] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (182 / 2) + 233] + fd_edgeFaceDst[ctr_1 - (156 / 2) + 2*(272 / 2) + 204];
            }
            {
              fd_edgeFaceDst[-(156 / 2) + 207] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(156 / 2) + 220] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(156 / 2) + 219] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(132 / 2) + 202] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_edgeFaceDst[-(156 / 2) + 207];
              fd_edgeFaceDst[-(156 / 2) + 2*(272 / 2) + 207] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(156 / 2) + 219] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(156 / 2) + 220] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_edgeFaceDst[-(156 / 2) + 2*(272 / 2) + 207];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(182 / 2) + 221] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 234] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 217] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_edgeFaceDst[-(182 / 2) + 221];
            fd_edgeFaceDst[-(182 / 2) + (272 / 2) + 221] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 234] + fd_edgeFaceDst[-(182 / 2) + (272 / 2) + 221];
          }
          {
            fd_edgeFaceDst[-(182 / 2) + 222] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 218] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_edgeFaceDst[-(182 / 2) + 222];
            fd_edgeFaceDst[-(182 / 2) + (272 / 2) + 222] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_edgeFaceDst[-(182 / 2) + (272 / 2) + 222];
            fd_edgeFaceDst[-(182 / 2) + 2*(272 / 2) + 222] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_edgeFaceDst[-(182 / 2) + 2*(272 / 2) + 222];
          }
          {
            fd_edgeFaceDst[-(182 / 2) + 223] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(156 / 2) + 219] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_edgeFaceDst[-(182 / 2) + 223];
            fd_edgeFaceDst[-(182 / 2) + 2*(272 / 2) + 223] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(182 / 2) + 237] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_edgeFaceDst[-(182 / 2) + 2*(272 / 2) + 223];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(210 / 2) + 238] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 235] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 270] + fd_edgeFaceDst[-(210 / 2) + 238];
          fd_edgeFaceDst[-(210 / 2) + (272 / 2) + 238] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(240 / 2) + 271] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(240 / 2) + 270] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(210 / 2) + 252] + fd_edgeFaceDst[-(210 / 2) + (272 / 2) + 238];
        }
        {
          fd_edgeFaceDst[-(210 / 2) + 239] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(182 / 2) + 236] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 271] + fd_edgeFaceDst[-(210 / 2) + 239];
          fd_edgeFaceDst[-(210 / 2) + 2*(272 / 2) + 239] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(210 / 2) + 253] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(210 / 2) + 254] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(240 / 2) + 271] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(240 / 2) + 270] + fd_edgeFaceDst[-(210 / 2) + 2*(272 / 2) + 239];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (240 / 2) + 255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (240 / 2) + 271] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (240 / 2) + 270] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (210 / 2) + 253] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (272 / 2) + 288] + fd_edgeFaceDst[ctr_1 - (240 / 2) + 255];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_5(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (1056 / 2)];
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (1056 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1056 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1056 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(1056 / 2) + 31] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 31] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 32] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 65] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 64] + fd_edgeFaceDst[-(0 / 2) + 2*(1056 / 2) + 31];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 33] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 35] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 34] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 68] + fd_edgeFaceDst[-(2 / 2) + 33];
          fd_edgeFaceDst[-(2 / 2) + (1056 / 2) + 33] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 69] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 68] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 34] + fd_edgeFaceDst[-(2 / 2) + (1056 / 2) + 33];
        }
        for (int ctr_1 = 1; ctr_1 < 30; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 33] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 33];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (1056 / 2) + 33] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 69] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (1056 / 2) + 33];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 34] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 67] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 65] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 64] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 31] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 98] + fd_edgeFaceDst[-(2 / 2) + 63];
          fd_edgeFaceDst[-(2 / 2) + 2*(1056 / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 64] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 65] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 98] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 97] + fd_edgeFaceDst[-(2 / 2) + 2*(1056 / 2) + 63];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 66] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 69] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 68] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 35] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 102] + fd_edgeFaceDst[-(6 / 2) + 66];
            fd_edgeFaceDst[-(6 / 2) + (1056 / 2) + 66] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 69] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 103] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 102] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 68] + fd_edgeFaceDst[-(6 / 2) + (1056 / 2) + 66];
          }
          for (int ctr_1 = 1; ctr_1 < 29; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 66] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 69] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 35] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 66];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (1056 / 2) + 66] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 69] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 103] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_edgeFaceDst[ctr_1 - (6 / 2) + (1056 / 2) + 66];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1056 / 2) + 66] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 68] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 69] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 101] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1056 / 2) + 66];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 95] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 98] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 97] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 64] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 131] + fd_edgeFaceDst[-(6 / 2) + 95];
            fd_edgeFaceDst[-(6 / 2) + 2*(1056 / 2) + 95] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 97] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 98] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 131] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 130] + fd_edgeFaceDst[-(6 / 2) + 2*(1056 / 2) + 95];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 99] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 103] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 102] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 69] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 136] + fd_edgeFaceDst[-(12 / 2) + 99];
              fd_edgeFaceDst[-(12 / 2) + (1056 / 2) + 99] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 103] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 137] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 136] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 102] + fd_edgeFaceDst[-(12 / 2) + (1056 / 2) + 99];
            }
            for (int ctr_1 = 1; ctr_1 < 28; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 99] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 103] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 69] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 99];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (1056 / 2) + 99] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 103] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 137] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_edgeFaceDst[ctr_1 - (12 / 2) + (1056 / 2) + 99];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(1056 / 2) + 99] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 102] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 103] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 135] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(1056 / 2) + 99];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 131] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 130] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 97] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 164] + fd_edgeFaceDst[-(12 / 2) + 127];
              fd_edgeFaceDst[-(12 / 2) + 2*(1056 / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 131] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 164] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 163] + fd_edgeFaceDst[-(12 / 2) + 2*(1056 / 2) + 127];
            }
          }
          {
            for (int ctr_1 = 0; ctr_1 < 28; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (20 / 2) + 132] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 137] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 103] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 170] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 132];
              if (ctr_1 + 4 < 31)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + (1056 / 2) + 132] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 137] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 171] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 170] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_edgeFaceDst[ctr_1 - (20 / 2) + (1056 / 2) + 132];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(1056 / 2) + 132] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 136] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 137] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 170] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 169] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(1056 / 2) + 132];
              }
            }
            for (int ctr_2 = 5; ctr_2 < 27; ctr_2 += 1)
              for (int ctr_1 = 0; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
                if (ctr_1 + ctr_2 < 31)
                {
                  fd_edgeFaceDst[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 35] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                }
                if (ctr_1 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                }
              }
            for (int ctr_1 = 0; ctr_1 < 5; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (756 / 2) + 891] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 919] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 918] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (702 / 2) + 885] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_edgeFaceDst[ctr_1 - (756 / 2) + 891];
              if (ctr_1 + 27 < 31)
              {
                fd_edgeFaceDst[ctr_1 - (756 / 2) + (1056 / 2) + 891] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 919] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 953] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 918] + fd_edgeFaceDst[ctr_1 - (756 / 2) + (1056 / 2) + 891];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 - (756 / 2) + 2*(1056 / 2) + 891] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 918] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 919] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 951] + fd_edgeFaceDst[ctr_1 - (756 / 2) + 2*(1056 / 2) + 891];
              }
            }
          }
          {
            {
              fd_edgeFaceDst[-(812 / 2) + 924] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(812 / 2) + 953] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(812 / 2) + 952] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(756 / 2) + 919] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(870 / 2) + 986] + fd_edgeFaceDst[-(812 / 2) + 924];
              fd_edgeFaceDst[-(812 / 2) + (1056 / 2) + 924] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(812 / 2) + 953] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(870 / 2) + 986] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(812 / 2) + 952] + fd_edgeFaceDst[-(812 / 2) + (1056 / 2) + 924];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (812 / 2) + 924] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 953] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (756 / 2) + 919] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (870 / 2) + 986] + fd_edgeFaceDst[ctr_1 - (812 / 2) + 924];
              fd_edgeFaceDst[ctr_1 - (812 / 2) + (1056 / 2) + 924] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 953] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (870 / 2) + 987] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (870 / 2) + 986] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_edgeFaceDst[ctr_1 - (812 / 2) + (1056 / 2) + 924];
              fd_edgeFaceDst[ctr_1 - (812 / 2) + 2*(1056 / 2) + 924] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 952] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (812 / 2) + 953] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (870 / 2) + 986] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (870 / 2) + 985] + fd_edgeFaceDst[ctr_1 - (812 / 2) + 2*(1056 / 2) + 924];
            }
            {
              fd_edgeFaceDst[-(812 / 2) + 927] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(812 / 2) + 956] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(812 / 2) + 955] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(756 / 2) + 922] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(870 / 2) + 989] + fd_edgeFaceDst[-(812 / 2) + 927];
              fd_edgeFaceDst[-(812 / 2) + 2*(1056 / 2) + 927] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(812 / 2) + 955] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(812 / 2) + 956] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(870 / 2) + 989] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_edgeFaceDst[-(812 / 2) + 2*(1056 / 2) + 927];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(870 / 2) + 957] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 986] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(812 / 2) + 953] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1020] + fd_edgeFaceDst[-(870 / 2) + 957];
            fd_edgeFaceDst[-(870 / 2) + (1056 / 2) + 957] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1020] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(870 / 2) + 986] + fd_edgeFaceDst[-(870 / 2) + (1056 / 2) + 957];
          }
          {
            fd_edgeFaceDst[-(870 / 2) + 958] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(812 / 2) + 954] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_edgeFaceDst[-(870 / 2) + 958];
            fd_edgeFaceDst[-(870 / 2) + (1056 / 2) + 958] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(930 / 2) + 1022] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_edgeFaceDst[-(870 / 2) + (1056 / 2) + 958];
            fd_edgeFaceDst[-(870 / 2) + 2*(1056 / 2) + 958] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1020] + fd_edgeFaceDst[-(870 / 2) + 2*(1056 / 2) + 958];
          }
          {
            fd_edgeFaceDst[-(870 / 2) + 959] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 989] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(812 / 2) + 955] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1022] + fd_edgeFaceDst[-(870 / 2) + 959];
            fd_edgeFaceDst[-(870 / 2) + 2*(1056 / 2) + 959] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(870 / 2) + 989] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(930 / 2) + 1022] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_edgeFaceDst[-(870 / 2) + 2*(1056 / 2) + 959];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(930 / 2) + 990] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(930 / 2) + 1020] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(870 / 2) + 987] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(992 / 2) + 1054] + fd_edgeFaceDst[-(930 / 2) + 990];
          fd_edgeFaceDst[-(930 / 2) + (1056 / 2) + 990] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(992 / 2) + 1055] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(992 / 2) + 1054] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(930 / 2) + 1020] + fd_edgeFaceDst[-(930 / 2) + (1056 / 2) + 990];
        }
        {
          fd_edgeFaceDst[-(930 / 2) + 991] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(930 / 2) + 1022] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(870 / 2) + 988] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(992 / 2) + 1055] + fd_edgeFaceDst[-(930 / 2) + 991];
          fd_edgeFaceDst[-(930 / 2) + 2*(1056 / 2) + 991] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(930 / 2) + 1021] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(930 / 2) + 1022] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(992 / 2) + 1055] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(992 / 2) + 1054] + fd_edgeFaceDst[-(930 / 2) + 2*(1056 / 2) + 991];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (992 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (992 / 2) + 1055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (992 / 2) + 1054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (930 / 2) + 1021] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1056 / 2) + 1088] + fd_edgeFaceDst[ctr_1 - (992 / 2) + 1023];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_6(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4160 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (4160 / 2)];
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4160 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (4160 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4160 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 65] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4160 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(4160 / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 63] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 64] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_edgeFaceDst[-(0 / 2) + 2*(4160 / 2) + 63];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 65] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 132] + fd_edgeFaceDst[-(2 / 2) + 65];
          fd_edgeFaceDst[-(2 / 2) + (4160 / 2) + 65] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 132] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 66] + fd_edgeFaceDst[-(2 / 2) + (4160 / 2) + 65];
        }
        for (int ctr_1 = 1; ctr_1 < 62; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 65] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 65];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (4160 / 2) + 65] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (4160 / 2) + 65];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(4160 / 2) + 65] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 131] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(4160 / 2) + 65];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 63] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 194] + fd_edgeFaceDst[-(2 / 2) + 127];
          fd_edgeFaceDst[-(2 / 2) + 2*(4160 / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 193] + fd_edgeFaceDst[-(2 / 2) + 2*(4160 / 2) + 127];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 130] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 133] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 132] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 198] + fd_edgeFaceDst[-(6 / 2) + 130];
            fd_edgeFaceDst[-(6 / 2) + (4160 / 2) + 130] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 199] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 198] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 132] + fd_edgeFaceDst[-(6 / 2) + (4160 / 2) + 130];
          }
          for (int ctr_1 = 1; ctr_1 < 61; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 130] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 133] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 67] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 198] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 130];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (4160 / 2) + 130] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 133] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 199] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 198] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_edgeFaceDst[ctr_1 - (6 / 2) + (4160 / 2) + 130];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(4160 / 2) + 130] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 132] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 133] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 198] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 197] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(4160 / 2) + 130];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 194] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 193] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 128] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 259] + fd_edgeFaceDst[-(6 / 2) + 191];
            fd_edgeFaceDst[-(6 / 2) + 2*(4160 / 2) + 191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 193] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 194] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 259] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 258] + fd_edgeFaceDst[-(6 / 2) + 2*(4160 / 2) + 191];
          }
        }
        for (int ctr_2 = 3; ctr_2 < 61; ctr_2 += 1)
        {
          {
            fd_edgeFaceDst[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_edgeFaceDst[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          }
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 67] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          }
          {
            fd_edgeFaceDst[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 64] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_edgeFaceDst[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63];
            fd_edgeFaceDst[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 64] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128] + fd_edgeFaceDst[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63];
          }
        }
        {
          {
            fd_edgeFaceDst[-(3782 / 2) + 3965] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3782 / 2) + 4027] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3782 / 2) + 4026] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3660 / 2) + 3961] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(3906 / 2) + 4092] + fd_edgeFaceDst[-(3782 / 2) + 3965];
            fd_edgeFaceDst[-(3782 / 2) + (4160 / 2) + 3965] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(3782 / 2) + 4027] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(3906 / 2) + 4092] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(3782 / 2) + 4026] + fd_edgeFaceDst[-(3782 / 2) + (4160 / 2) + 3965];
          }
          for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (3782 / 2) + 3965] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4027] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4026] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3660 / 2) + 3961] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4092] + fd_edgeFaceDst[ctr_1 - (3782 / 2) + 3965];
            fd_edgeFaceDst[ctr_1 - (3782 / 2) + (4160 / 2) + 3965] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4027] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4093] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4092] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4026] + fd_edgeFaceDst[ctr_1 - (3782 / 2) + (4160 / 2) + 3965];
            fd_edgeFaceDst[ctr_1 - (3782 / 2) + 2*(4160 / 2) + 3965] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4026] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (3782 / 2) + 4027] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4092] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4091] + fd_edgeFaceDst[ctr_1 - (3782 / 2) + 2*(4160 / 2) + 3965];
          }
          {
            fd_edgeFaceDst[-(3782 / 2) + 3967] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3782 / 2) + 4029] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3782 / 2) + 4028] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3660 / 2) + 3963] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_edgeFaceDst[-(3782 / 2) + 3967];
            fd_edgeFaceDst[-(3782 / 2) + 2*(4160 / 2) + 3967] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(3782 / 2) + 4028] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(3782 / 2) + 4029] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_edgeFaceDst[-(3782 / 2) + 2*(4160 / 2) + 3967];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(3906 / 2) + 4030] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4092] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3782 / 2) + 4027] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4158] + fd_edgeFaceDst[-(3906 / 2) + 4030];
          fd_edgeFaceDst[-(3906 / 2) + (4160 / 2) + 4030] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(4032 / 2) + 4159] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(4032 / 2) + 4158] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(3906 / 2) + 4092] + fd_edgeFaceDst[-(3906 / 2) + (4160 / 2) + 4030];
        }
        {
          fd_edgeFaceDst[-(3906 / 2) + 4031] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(3782 / 2) + 4028] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4159] + fd_edgeFaceDst[-(3906 / 2) + 4031];
          fd_edgeFaceDst[-(3906 / 2) + 2*(4160 / 2) + 4031] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(3906 / 2) + 4093] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(3906 / 2) + 4094] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(4032 / 2) + 4159] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(4032 / 2) + 4158] + fd_edgeFaceDst[-(3906 / 2) + 2*(4160 / 2) + 4031];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (4032 / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4032 / 2) + 4159] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4032 / 2) + 4158] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (3906 / 2) + 4093] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4160 / 2) + 4224] + fd_edgeFaceDst[ctr_1 - (4032 / 2) + 4095];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_7(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16512 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (16512 / 2)];
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16512 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (16512 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16512 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 129] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16512 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(16512 / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 127] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 257] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 256] + fd_edgeFaceDst[-(0 / 2) + 2*(16512 / 2) + 127];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 129] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 131] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 130] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 260] + fd_edgeFaceDst[-(2 / 2) + 129];
          fd_edgeFaceDst[-(2 / 2) + (16512 / 2) + 129] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 261] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 260] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 130] + fd_edgeFaceDst[-(2 / 2) + (16512 / 2) + 129];
        }
        for (int ctr_1 = 1; ctr_1 < 126; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 129] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 129];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (16512 / 2) + 129] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 261] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (16512 / 2) + 129];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(16512 / 2) + 129] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 259] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(16512 / 2) + 129];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 257] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 256] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 127] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 386] + fd_edgeFaceDst[-(2 / 2) + 255];
          fd_edgeFaceDst[-(2 / 2) + 2*(16512 / 2) + 255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 256] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 257] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 386] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 385] + fd_edgeFaceDst[-(2 / 2) + 2*(16512 / 2) + 255];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 258] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 261] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 260] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 131] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 390] + fd_edgeFaceDst[-(6 / 2) + 258];
            fd_edgeFaceDst[-(6 / 2) + (16512 / 2) + 258] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 261] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 391] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 390] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 260] + fd_edgeFaceDst[-(6 / 2) + (16512 / 2) + 258];
          }
          for (int ctr_1 = 1; ctr_1 < 125; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 258] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 261] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 131] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 258];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (16512 / 2) + 258] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 261] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 391] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_edgeFaceDst[ctr_1 - (6 / 2) + (16512 / 2) + 258];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(16512 / 2) + 258] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 260] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 261] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 389] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(16512 / 2) + 258];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 386] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 385] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 256] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 515] + fd_edgeFaceDst[-(6 / 2) + 383];
            fd_edgeFaceDst[-(6 / 2) + 2*(16512 / 2) + 383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 385] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 386] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 515] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 514] + fd_edgeFaceDst[-(6 / 2) + 2*(16512 / 2) + 383];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 387] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 391] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 390] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 261] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 520] + fd_edgeFaceDst[-(12 / 2) + 387];
              fd_edgeFaceDst[-(12 / 2) + (16512 / 2) + 387] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 391] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 521] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 520] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 390] + fd_edgeFaceDst[-(12 / 2) + (16512 / 2) + 387];
            }
            for (int ctr_1 = 1; ctr_1 < 124; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 387] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 391] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 261] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 520] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 387];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (16512 / 2) + 387] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 391] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 521] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 520] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_edgeFaceDst[ctr_1 - (12 / 2) + (16512 / 2) + 387];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(16512 / 2) + 387] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 390] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 391] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 520] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 519] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(16512 / 2) + 387];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 511] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 515] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 514] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 385] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 644] + fd_edgeFaceDst[-(12 / 2) + 511];
              fd_edgeFaceDst[-(12 / 2) + 2*(16512 / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 515] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 644] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 643] + fd_edgeFaceDst[-(12 / 2) + 2*(16512 / 2) + 511];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 124; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_edgeFaceDst[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 131] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            }
            {
              fd_edgeFaceDst[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 128] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_edgeFaceDst[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127];
              fd_edgeFaceDst[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 128] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256] + fd_edgeFaceDst[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127];
            }
          }
          {
            {
              fd_edgeFaceDst[-(15500 / 2) + 15996] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(15500 / 2) + 16121] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(15500 / 2) + 16120] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15252 / 2) + 15991] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(15750 / 2) + 16250] + fd_edgeFaceDst[-(15500 / 2) + 15996];
              fd_edgeFaceDst[-(15500 / 2) + (16512 / 2) + 15996] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(15500 / 2) + 16121] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(15750 / 2) + 16250] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(15500 / 2) + 16120] + fd_edgeFaceDst[-(15500 / 2) + (16512 / 2) + 15996];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (15500 / 2) + 15996] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16121] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16120] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (15252 / 2) + 15991] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (15750 / 2) + 16250] + fd_edgeFaceDst[ctr_1 - (15500 / 2) + 15996];
              fd_edgeFaceDst[ctr_1 - (15500 / 2) + (16512 / 2) + 15996] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16121] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (15750 / 2) + 16251] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (15750 / 2) + 16250] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16120] + fd_edgeFaceDst[ctr_1 - (15500 / 2) + (16512 / 2) + 15996];
              fd_edgeFaceDst[ctr_1 - (15500 / 2) + 2*(16512 / 2) + 15996] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16120] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (15500 / 2) + 16121] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (15750 / 2) + 16250] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (15750 / 2) + 16249] + fd_edgeFaceDst[ctr_1 - (15500 / 2) + 2*(16512 / 2) + 15996];
            }
            {
              fd_edgeFaceDst[-(15500 / 2) + 15999] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(15500 / 2) + 16124] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(15500 / 2) + 16123] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15252 / 2) + 15994] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(15750 / 2) + 16253] + fd_edgeFaceDst[-(15500 / 2) + 15999];
              fd_edgeFaceDst[-(15500 / 2) + 2*(16512 / 2) + 15999] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(15500 / 2) + 16123] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(15500 / 2) + 16124] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(15750 / 2) + 16253] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_edgeFaceDst[-(15500 / 2) + 2*(16512 / 2) + 15999];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(15750 / 2) + 16125] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16250] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15500 / 2) + 16121] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16380] + fd_edgeFaceDst[-(15750 / 2) + 16125];
            fd_edgeFaceDst[-(15750 / 2) + (16512 / 2) + 16125] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16380] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(15750 / 2) + 16250] + fd_edgeFaceDst[-(15750 / 2) + (16512 / 2) + 16125];
          }
          {
            fd_edgeFaceDst[-(15750 / 2) + 16126] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15500 / 2) + 16122] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_edgeFaceDst[-(15750 / 2) + 16126];
            fd_edgeFaceDst[-(15750 / 2) + (16512 / 2) + 16126] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16002 / 2) + 16382] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_edgeFaceDst[-(15750 / 2) + (16512 / 2) + 16126];
            fd_edgeFaceDst[-(15750 / 2) + 2*(16512 / 2) + 16126] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16380] + fd_edgeFaceDst[-(15750 / 2) + 2*(16512 / 2) + 16126];
          }
          {
            fd_edgeFaceDst[-(15750 / 2) + 16127] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16253] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15500 / 2) + 16123] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16382] + fd_edgeFaceDst[-(15750 / 2) + 16127];
            fd_edgeFaceDst[-(15750 / 2) + 2*(16512 / 2) + 16127] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(15750 / 2) + 16253] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16002 / 2) + 16382] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_edgeFaceDst[-(15750 / 2) + 2*(16512 / 2) + 16127];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(16002 / 2) + 16254] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16002 / 2) + 16380] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15750 / 2) + 16251] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16256 / 2) + 16510] + fd_edgeFaceDst[-(16002 / 2) + 16254];
          fd_edgeFaceDst[-(16002 / 2) + (16512 / 2) + 16254] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16256 / 2) + 16511] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16256 / 2) + 16510] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16002 / 2) + 16380] + fd_edgeFaceDst[-(16002 / 2) + (16512 / 2) + 16254];
        }
        {
          fd_edgeFaceDst[-(16002 / 2) + 16255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16002 / 2) + 16382] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(15750 / 2) + 16252] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16256 / 2) + 16511] + fd_edgeFaceDst[-(16002 / 2) + 16255];
          fd_edgeFaceDst[-(16002 / 2) + 2*(16512 / 2) + 16255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16002 / 2) + 16381] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16002 / 2) + 16382] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16256 / 2) + 16511] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16256 / 2) + 16510] + fd_edgeFaceDst[-(16002 / 2) + 2*(16512 / 2) + 16255];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (16256 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16256 / 2) + 16511] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16256 / 2) + 16510] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16002 / 2) + 16381] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16512 / 2) + 16640] + fd_edgeFaceDst[ctr_1 - (16256 / 2) + 16383];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_8(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (65792 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (65792 / 2)];
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (65792 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (65792 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(65792 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 258] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 257] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(65792 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(65792 / 2) + 255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 255] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 256] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 513] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 512] + fd_edgeFaceDst[-(0 / 2) + 2*(65792 / 2) + 255];
    }
    for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_edgeFaceDst[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 259] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 256] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_edgeFaceDst[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255];
        fd_edgeFaceDst[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 256] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512] + fd_edgeFaceDst[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255];
      }
    }
    {
      fd_edgeFaceDst[-(65280 / 2) + 65535] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(65280 / 2) + 65791] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(65280 / 2) + 65790] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(65792 / 2) + 66048] + fd_edgeFaceDst[-(65280 / 2) + 65535];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (65280 / 2) + 65535] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (65280 / 2) + 65791] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (65280 / 2) + 65790] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (65792 / 2) + 66048] + fd_edgeFaceDst[ctr_1 - (65280 / 2) + 65535];
      }
      fd_edgeFaceDst[-(65280 / 2) + 65535] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(65280 / 2) + 65791] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(65280 / 2) + 65790] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(64770 / 2) + 65533] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(65792 / 2) + 66048] + fd_edgeFaceDst[-(65280 / 2) + 65535];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_9(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (262656 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (262656 / 2)];
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (262656 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (262656 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(262656 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(262656 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(262656 / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 511] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 512] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 1025] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 1024] + fd_edgeFaceDst[-(0 / 2) + 2*(262656 / 2) + 511];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 513] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 515] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 514] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 1028] + fd_edgeFaceDst[-(2 / 2) + 513];
          fd_edgeFaceDst[-(2 / 2) + (262656 / 2) + 513] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 1029] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 1028] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 514] + fd_edgeFaceDst[-(2 / 2) + (262656 / 2) + 513];
        }
        for (int ctr_1 = 1; ctr_1 < 510; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 513] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 515] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 1028] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 513];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (262656 / 2) + 513] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 1029] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 1028] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (262656 / 2) + 513];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 515] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 1028] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 1027] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 1025] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 1024] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 511] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 1538] + fd_edgeFaceDst[-(2 / 2) + 1023];
          fd_edgeFaceDst[-(2 / 2) + 2*(262656 / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 1024] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 1025] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 1538] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 1537] + fd_edgeFaceDst[-(2 / 2) + 2*(262656 / 2) + 1023];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 510; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_edgeFaceDst[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 515] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        }
        {
          fd_edgeFaceDst[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 512] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_edgeFaceDst[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511];
          fd_edgeFaceDst[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 512] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024] + fd_edgeFaceDst[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511];
        }
      }
      {
        {
          fd_edgeFaceDst[-(260610 / 2) + 261630] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(260610 / 2) + 262140] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(259590 / 2) + 261627] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(261632 / 2) + 262654] + fd_edgeFaceDst[-(260610 / 2) + 261630];
          fd_edgeFaceDst[-(260610 / 2) + (262656 / 2) + 261630] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(261632 / 2) + 262655] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(261632 / 2) + 262654] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(260610 / 2) + 262140] + fd_edgeFaceDst[-(260610 / 2) + (262656 / 2) + 261630];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (260610 / 2) + 261630] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (260610 / 2) + 262140] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (259590 / 2) + 261627] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (261632 / 2) + 262654] + fd_edgeFaceDst[ctr_1 - (260610 / 2) + 261630];
        }
        {
          fd_edgeFaceDst[-(260610 / 2) + 261631] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(260610 / 2) + 262142] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(259590 / 2) + 261628] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(261632 / 2) + 262655] + fd_edgeFaceDst[-(260610 / 2) + 261631];
          fd_edgeFaceDst[-(260610 / 2) + 2*(262656 / 2) + 261631] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(260610 / 2) + 262141] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(260610 / 2) + 262142] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(261632 / 2) + 262655] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(261632 / 2) + 262654] + fd_edgeFaceDst[-(260610 / 2) + 2*(262656 / 2) + 261631];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (261632 / 2) + 262143] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (261632 / 2) + 262655] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (261632 / 2) + 262654] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (260610 / 2) + 262141] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (262656 / 2) + 263168] + fd_edgeFaceDst[ctr_1 - (261632 / 2) + 262143];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_10(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1049600 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (1049600 / 2)];
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1049600 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (1049600 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1049600 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1049600 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(1049600 / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1023] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 1024] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2049] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 2048] + fd_edgeFaceDst[-(0 / 2) + 2*(1049600 / 2) + 1023];
    }
    for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_edgeFaceDst[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1027] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1024] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_edgeFaceDst[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023];
        fd_edgeFaceDst[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1024] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2048] + fd_edgeFaceDst[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023];
      }
    }
    {
      fd_edgeFaceDst[-(1047552 / 2) + 1048575] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1047552 / 2) + 1049599] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1047552 / 2) + 1049598] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1049600 / 2) + 1050624] + fd_edgeFaceDst[-(1047552 / 2) + 1048575];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (1047552 / 2) + 1048575] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (1047552 / 2) + 1049599] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (1047552 / 2) + 1049598] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (1049600 / 2) + 1050624] + fd_edgeFaceDst[ctr_1 - (1047552 / 2) + 1048575];
      }
      fd_edgeFaceDst[-(1047552 / 2) + 1048575] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(1047552 / 2) + 1049599] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(1047552 / 2) + 1049598] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(1045506 / 2) + 1048573] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(1049600 / 2) + 1050624] + fd_edgeFaceDst[-(1047552 / 2) + 1048575];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_11(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4196352 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (4196352 / 2)];
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4196352 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (4196352 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4196352 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2050] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 2049] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4196352 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(4196352 / 2) + 2047] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 2047] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 2048] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4097] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 4096] + fd_edgeFaceDst[-(0 / 2) + 2*(4196352 / 2) + 2047];
    }
    for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_edgeFaceDst[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2051] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2048] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_edgeFaceDst[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047];
        fd_edgeFaceDst[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2048] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096] + fd_edgeFaceDst[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047];
      }
    }
    {
      fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4192256 / 2) + 4196351] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4192256 / 2) + 4196350] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4196352 / 2) + 4198400] + fd_edgeFaceDst[-(4192256 / 2) + 4194303];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (4192256 / 2) + 4194303] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (4192256 / 2) + 4196351] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (4192256 / 2) + 4196350] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (4196352 / 2) + 4198400] + fd_edgeFaceDst[ctr_1 - (4192256 / 2) + 4194303];
      }
      fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(4192256 / 2) + 4196351] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(4192256 / 2) + 4196350] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(4188162 / 2) + 4194301] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(4196352 / 2) + 4198400] + fd_edgeFaceDst[-(4192256 / 2) + 4194303];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_12(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16781312 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (16781312 / 2)];
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16781312 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (16781312 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16781312 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16781312 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(16781312 / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 4095] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 4096] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8193] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 8192] + fd_edgeFaceDst[-(0 / 2) + 2*(16781312 / 2) + 4095];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 4097] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 4099] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 4098] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 8196] + fd_edgeFaceDst[-(2 / 2) + 4097];
          fd_edgeFaceDst[-(2 / 2) + (16781312 / 2) + 4097] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 8197] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 8196] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 4098] + fd_edgeFaceDst[-(2 / 2) + (16781312 / 2) + 4097];
        }
        for (int ctr_1 = 1; ctr_1 < 4094; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 4097] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 4097];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8197] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_edgeFaceDst[ctr_1 - (2 / 2) + (16781312 / 2) + 4097];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8195] + fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 8191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8193] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8192] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(0 / 2) + 4095] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 12290] + fd_edgeFaceDst[-(2 / 2) + 8191];
          fd_edgeFaceDst[-(2 / 2) + 2*(16781312 / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(2 / 2) + 8192] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8193] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 12290] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 12289] + fd_edgeFaceDst[-(2 / 2) + 2*(16781312 / 2) + 8191];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 8194] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 8197] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 8196] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 4099] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 12294] + fd_edgeFaceDst[-(6 / 2) + 8194];
            fd_edgeFaceDst[-(6 / 2) + (16781312 / 2) + 8194] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 8197] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 12295] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 12294] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(6 / 2) + 8196] + fd_edgeFaceDst[-(6 / 2) + (16781312 / 2) + 8194];
          }
          for (int ctr_1 = 1; ctr_1 < 4093; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 8194] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8197] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 4099] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 8194];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (16781312 / 2) + 8194] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8197] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12295] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_edgeFaceDst[ctr_1 - (6 / 2) + (16781312 / 2) + 8194];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(16781312 / 2) + 8194] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8196] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8197] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12293] + fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(16781312 / 2) + 8194];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 12287] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 12290] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 12289] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8192] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 16387] + fd_edgeFaceDst[-(6 / 2) + 12287];
            fd_edgeFaceDst[-(6 / 2) + 2*(16781312 / 2) + 12287] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(6 / 2) + 12289] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(6 / 2) + 12290] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 16387] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 16386] + fd_edgeFaceDst[-(6 / 2) + 2*(16781312 / 2) + 12287];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 12291] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 12295] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 12294] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 8197] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 16392] + fd_edgeFaceDst[-(12 / 2) + 12291];
              fd_edgeFaceDst[-(12 / 2) + (16781312 / 2) + 12291] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 12295] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 16393] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 16392] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(12 / 2) + 12294] + fd_edgeFaceDst[-(12 / 2) + (16781312 / 2) + 12291];
            }
            for (int ctr_1 = 1; ctr_1 < 4092; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 12291] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12295] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (6 / 2) + 8197] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 12291];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (16781312 / 2) + 12291] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12295] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16393] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_edgeFaceDst[ctr_1 - (12 / 2) + (16781312 / 2) + 12291];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(16781312 / 2) + 12291] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12294] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12295] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16391] + fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(16781312 / 2) + 12291];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 16387] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 16386] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(6 / 2) + 12289] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 20484] + fd_edgeFaceDst[-(12 / 2) + 16383];
              fd_edgeFaceDst[-(12 / 2) + 2*(16781312 / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(12 / 2) + 16386] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(12 / 2) + 16387] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(20 / 2) + 20484] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 20483] + fd_edgeFaceDst[-(12 / 2) + 2*(16781312 / 2) + 16383];
            }
          }
          {
            {
              {
                fd_edgeFaceDst[-(20 / 2) + 16388] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 16393] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 16392] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 12295] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 20490] + fd_edgeFaceDst[-(20 / 2) + 16388];
                fd_edgeFaceDst[-(20 / 2) + (16781312 / 2) + 16388] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 16393] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(30 / 2) + 20491] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 20490] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(20 / 2) + 16392] + fd_edgeFaceDst[-(20 / 2) + (16781312 / 2) + 16388];
              }
              for (int ctr_1 = 1; ctr_1 < 4091; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 16388] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16393] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (12 / 2) + 12295] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 20490] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 16388];
                fd_edgeFaceDst[ctr_1 - (20 / 2) + (16781312 / 2) + 16388] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16393] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 20491] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 20490] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_edgeFaceDst[ctr_1 - (20 / 2) + (16781312 / 2) + 16388];
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(16781312 / 2) + 16388] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16392] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (20 / 2) + 16393] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 20490] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (30 / 2) + 20489] + fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(16781312 / 2) + 16388];
              }
              {
                fd_edgeFaceDst[-(20 / 2) + 20479] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 20484] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 20483] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(12 / 2) + 16386] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 24581] + fd_edgeFaceDst[-(20 / 2) + 20479];
                fd_edgeFaceDst[-(20 / 2) + 2*(16781312 / 2) + 20479] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(20 / 2) + 20483] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(20 / 2) + 20484] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(30 / 2) + 24581] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(30 / 2) + 24580] + fd_edgeFaceDst[-(20 / 2) + 2*(16781312 / 2) + 20479];
              }
            }
            for (int ctr_2 = 5; ctr_2 < 4091; ctr_2 += 1)
            {
              {
                fd_edgeFaceDst[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_edgeFaceDst[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              }
              for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4099] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              }
              {
                fd_edgeFaceDst[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4096] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_edgeFaceDst[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095];
                fd_edgeFaceDst[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4096] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192] + fd_edgeFaceDst[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095];
              }
            }
            {
              {
                fd_edgeFaceDst[-(16740372 / 2) + 16760827] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16740372 / 2) + 16764919] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16740372 / 2) + 16764918] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16732190 / 2) + 16760821] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16748556 / 2) + 16769016] + fd_edgeFaceDst[-(16740372 / 2) + 16760827];
                fd_edgeFaceDst[-(16740372 / 2) + (16781312 / 2) + 16760827] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16740372 / 2) + 16764919] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16748556 / 2) + 16769017] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16748556 / 2) + 16769016] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16740372 / 2) + 16764918] + fd_edgeFaceDst[-(16740372 / 2) + (16781312 / 2) + 16760827];
              }
              for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 - (16740372 / 2) + 16760827] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764919] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764918] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16732190 / 2) + 16760821] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_edgeFaceDst[ctr_1 - (16740372 / 2) + 16760827];
                fd_edgeFaceDst[ctr_1 - (16740372 / 2) + (16781312 / 2) + 16760827] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764919] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769017] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764918] + fd_edgeFaceDst[ctr_1 - (16740372 / 2) + (16781312 / 2) + 16760827];
                fd_edgeFaceDst[ctr_1 - (16740372 / 2) + 2*(16781312 / 2) + 16760827] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764918] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764919] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769015] + fd_edgeFaceDst[ctr_1 - (16740372 / 2) + 2*(16781312 / 2) + 16760827];
              }
              {
                fd_edgeFaceDst[-(16740372 / 2) + 16760831] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16740372 / 2) + 16764923] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16740372 / 2) + 16764922] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16732190 / 2) + 16760825] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16748556 / 2) + 16769020] + fd_edgeFaceDst[-(16740372 / 2) + 16760831];
                fd_edgeFaceDst[-(16740372 / 2) + 2*(16781312 / 2) + 16760831] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16740372 / 2) + 16764922] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16740372 / 2) + 16764923] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16748556 / 2) + 16769020] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16748556 / 2) + 16769019] + fd_edgeFaceDst[-(16740372 / 2) + 2*(16781312 / 2) + 16760831];
              }
            }
          }
          {
            {
              fd_edgeFaceDst[-(16748556 / 2) + 16764924] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16748556 / 2) + 16769017] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16748556 / 2) + 16769016] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16740372 / 2) + 16764919] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16756742 / 2) + 16773114] + fd_edgeFaceDst[-(16748556 / 2) + 16764924];
              fd_edgeFaceDst[-(16748556 / 2) + (16781312 / 2) + 16764924] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16748556 / 2) + 16769017] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16756742 / 2) + 16773114] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16748556 / 2) + 16769016] + fd_edgeFaceDst[-(16748556 / 2) + (16781312 / 2) + 16764924];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (16748556 / 2) + 16764924] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769017] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16740372 / 2) + 16764919] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16756742 / 2) + 16773114] + fd_edgeFaceDst[ctr_1 - (16748556 / 2) + 16764924];
              fd_edgeFaceDst[ctr_1 - (16748556 / 2) + (16781312 / 2) + 16764924] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769017] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16756742 / 2) + 16773115] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16756742 / 2) + 16773114] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_edgeFaceDst[ctr_1 - (16748556 / 2) + (16781312 / 2) + 16764924];
              fd_edgeFaceDst[ctr_1 - (16748556 / 2) + 2*(16781312 / 2) + 16764924] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769016] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16748556 / 2) + 16769017] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16756742 / 2) + 16773114] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16756742 / 2) + 16773113] + fd_edgeFaceDst[ctr_1 - (16748556 / 2) + 2*(16781312 / 2) + 16764924];
            }
            {
              fd_edgeFaceDst[-(16748556 / 2) + 16764927] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16748556 / 2) + 16769020] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16748556 / 2) + 16769019] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16740372 / 2) + 16764922] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16756742 / 2) + 16773117] + fd_edgeFaceDst[-(16748556 / 2) + 16764927];
              fd_edgeFaceDst[-(16748556 / 2) + 2*(16781312 / 2) + 16764927] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16748556 / 2) + 16769019] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16748556 / 2) + 16769020] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16756742 / 2) + 16773117] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_edgeFaceDst[-(16748556 / 2) + 2*(16781312 / 2) + 16764927];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(16756742 / 2) + 16769021] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773114] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16748556 / 2) + 16769017] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777212] + fd_edgeFaceDst[-(16756742 / 2) + 16769021];
            fd_edgeFaceDst[-(16756742 / 2) + (16781312 / 2) + 16769021] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777212] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16756742 / 2) + 16773114] + fd_edgeFaceDst[-(16756742 / 2) + (16781312 / 2) + 16769021];
          }
          {
            fd_edgeFaceDst[-(16756742 / 2) + 16769022] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16748556 / 2) + 16769018] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_edgeFaceDst[-(16756742 / 2) + 16769022];
            fd_edgeFaceDst[-(16756742 / 2) + (16781312 / 2) + 16769022] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16764930 / 2) + 16777214] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_edgeFaceDst[-(16756742 / 2) + (16781312 / 2) + 16769022];
            fd_edgeFaceDst[-(16756742 / 2) + 2*(16781312 / 2) + 16769022] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777212] + fd_edgeFaceDst[-(16756742 / 2) + 2*(16781312 / 2) + 16769022];
          }
          {
            fd_edgeFaceDst[-(16756742 / 2) + 16769023] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773117] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16748556 / 2) + 16769019] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777214] + fd_edgeFaceDst[-(16756742 / 2) + 16769023];
            fd_edgeFaceDst[-(16756742 / 2) + 2*(16781312 / 2) + 16769023] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16756742 / 2) + 16773117] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16764930 / 2) + 16777214] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_edgeFaceDst[-(16756742 / 2) + 2*(16781312 / 2) + 16769023];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(16764930 / 2) + 16773118] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16764930 / 2) + 16777212] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16756742 / 2) + 16773115] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16773120 / 2) + 16781310] + fd_edgeFaceDst[-(16764930 / 2) + 16773118];
          fd_edgeFaceDst[-(16764930 / 2) + (16781312 / 2) + 16773118] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(16773120 / 2) + 16781311] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(16773120 / 2) + 16781310] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(16764930 / 2) + 16777212] + fd_edgeFaceDst[-(16764930 / 2) + (16781312 / 2) + 16773118];
        }
        {
          fd_edgeFaceDst[-(16764930 / 2) + 16773119] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(16764930 / 2) + 16777214] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(16756742 / 2) + 16773116] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(16773120 / 2) + 16781311] + fd_edgeFaceDst[-(16764930 / 2) + 16773119];
          fd_edgeFaceDst[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(16764930 / 2) + 16777213] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(16764930 / 2) + 16777214] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(16773120 / 2) + 16781311] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(16773120 / 2) + 16781310] + fd_edgeFaceDst[-(16764930 / 2) + 2*(16781312 / 2) + 16773119];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (16773120 / 2) + 16777215] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (16773120 / 2) + 16781311] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (16773120 / 2) + 16781310] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (16764930 / 2) + 16777213] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (16781312 / 2) + 16785408] + fd_edgeFaceDst[ctr_1 - (16773120 / 2) + 16777215];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_13(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (67117056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (67117056 / 2)];
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (67117056 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (67117056 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(67117056 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(67117056 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(67117056 / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 8191] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 8192] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 16385] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 16384] + fd_edgeFaceDst[-(0 / 2) + 2*(67117056 / 2) + 8191];
    }
    for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_edgeFaceDst[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8195] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8192] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_edgeFaceDst[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191];
        fd_edgeFaceDst[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8192] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384] + fd_edgeFaceDst[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191];
      }
    }
    {
      fd_edgeFaceDst[-(67100672 / 2) + 67108863] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67100672 / 2) + 67117055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67100672 / 2) + 67117054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67117056 / 2) + 67125248] + fd_edgeFaceDst[-(67100672 / 2) + 67108863];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (67100672 / 2) + 67108863] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (67100672 / 2) + 67117055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (67100672 / 2) + 67117054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (67117056 / 2) + 67125248] + fd_edgeFaceDst[ctr_1 - (67100672 / 2) + 67108863];
      }
      fd_edgeFaceDst[-(67100672 / 2) + 67108863] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(67100672 / 2) + 67117055] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(67100672 / 2) + 67117054] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(67084290 / 2) + 67108861] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(67117056 / 2) + 67125248] + fd_edgeFaceDst[-(67100672 / 2) + 67108863];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_14(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (268451840 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[-(2 / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[-(0 / 2)] + fd_edgeFaceDst[-(0 / 2) + (268451840 / 2)];
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (268451840 / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceDst[ctr_1 - (0 / 2) + (268451840 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(268451840 / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (0 / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (0 / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16386] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(268451840 / 2)];
      }
      fd_edgeFaceDst[-(0 / 2) + 2*(268451840 / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[-(0 / 2) + 16383] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[-(0 / 2) + 16384] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[-(2 / 2) + 32769] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[-(2 / 2) + 32768] + fd_edgeFaceDst[-(0 / 2) + 2*(268451840 / 2) + 16383];
    }
    for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_edgeFaceDst[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16387] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      {
        fd_edgeFaceDst[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16384] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32769] + fd_edgeFaceDst[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383];
        fd_edgeFaceDst[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16384] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32769] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32768] + fd_edgeFaceDst[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383];
      }
    }
    {
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(268451840 / 2) + 268468224] + fd_edgeFaceDst[-(268419072 / 2) + 268435455];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 - (268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 - (268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 - (268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 - (268451840 / 2) + 268468224] + fd_edgeFaceDst[ctr_1 - (268419072 / 2) + 268435455];
      }
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[-(268419072 / 2) + 268451839] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[-(268419072 / 2) + 268451838] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[-(268386306 / 2) + 268435453] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[-(268451840 / 2) + 268468224] + fd_edgeFaceDst[-(268419072 / 2) + 268435455];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_any(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  const double fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  const double fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < (1 << level); ctr_2 += 1)
    for (int ctr_1 = 0; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      if (ctr_2 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      if (ctr_1 + ctr_2 < (1 << level) - 1)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1] + fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
      }
      if (ctr_1 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
    }
}




static void apply_2D_macroface_vertexdof_to_edgedof_add(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  switch( level )
  {
#if 0
    case 2:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_2(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 3:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_3(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 4:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_4(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 5:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_5(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 6:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_6(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 7:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_7(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 8:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_8(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 9:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_9(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 10:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_10(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 11:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_11(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 12:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_12(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 13:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_13(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
    case 14:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_14(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil);
      break;
#endif
    default:
      apply_2D_macroface_vertexdof_to_edgedof_add_level_any(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil, level);
      break;
  }
}




void applyFaceReplace( double*          fd_edgeFaceDst,
                       double*          fd_vertexFaceSrc,
                       double*          fd_vertexToEdgeFaceStencil,
                       walberla::uint_t level )
{
  apply_2D_macroface_vertexdof_to_edgedof_replace(fd_edgeFaceDst, fd_vertexFaceSrc,
      &fd_vertexToEdgeFaceStencil[4],
      &fd_vertexToEdgeFaceStencil[0],
      &fd_vertexToEdgeFaceStencil[8],
      static_cast< int64_t >(level));
}

void applyFaceAdd( double*          fd_edgeFaceDst,
                   double*          fd_vertexFaceSrc,
                   double*          fd_vertexToEdgeFaceStencil,
                   walberla::uint_t level )
{
  apply_2D_macroface_vertexdof_to_edgedof_add(fd_edgeFaceDst, fd_vertexFaceSrc,
      &fd_vertexToEdgeFaceStencil[4],
      &fd_vertexToEdgeFaceStencil[0],
      &fd_vertexToEdgeFaceStencil[8],
      static_cast< int64_t >(level));
}


}
}
}


#include "generatedKernels.hpp"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {


static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_2(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5];
    }
    for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5];
    }
    for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5];
    }
  }
  for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_3(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9];
    }
    for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9];
    }
    for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9];
    }
  }
  for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_4(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17];
    }
    for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17];
    }
    for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17];
    }
  }
  for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_5(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33];
    }
    for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33];
    }
    for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33];
    }
  }
  for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_6(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65];
    }
    for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65];
    }
    for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65];
    }
  }
  for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_7(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129];
    }
    for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129];
    }
    for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129];
    }
  }
  for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_8(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257];
    }
    for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257];
    }
    for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257];
    }
  }
  for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_9(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513];
    }
    for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513];
    }
    for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513];
    }
  }
  for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_10(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025];
    }
    for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025];
    }
    for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025];
    }
  }
  for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_11(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049];
    }
    for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049];
    }
    for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049];
    }
  }
  for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_12(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097];
    }
    for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097];
    }
    for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097];
    }
  }
  for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_13(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193];
    }
    for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193];
    }
    for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193];
    }
  }
  for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_14(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385];
    }
    for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385];
    }
    for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385];
    }
  }
  for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1];
    }
    for (int ctr_1 = (1 << level) - 1; ctr_1 < (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < (1 << level) - 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + (1 << level) - 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1];
    }
  }
  for (int ctr_2 = (1 << level) - 1; ctr_2 < (1 << level); ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)];
    }
  }
}




static void apply_2D_macroface_vertexdof_to_edgedof_replace(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level)
{
  switch( level )
  {
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
    default:
      apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(fd_edgeFaceDst, fd_vertexFaceSrc, fd_vertexToDiagonalEdgeFaceStencil, fd_vertexToHorizontalEdgeFaceStencil, fd_vertexToVerticalEdgeFaceStencil, level);
      break;
  }
}



/// Add start

static void apply_2D_macroface_vertexdof_to_edgedof_add_level_2(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 7] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 6] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_3(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 11] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 10] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_4(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 19] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 18] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_5(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 35] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 34] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_6(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 67] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 66] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_7(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 131] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 130] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_8(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 259] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 258] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_9(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 515] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 514] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_10(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1027] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1026] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_11(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2051] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2050] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_12(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4099] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4098] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_13(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8195] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8194] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_vertexdof_to_edgedof_add_level_14(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil)
{
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil0 = fd_vertexToDiagonalEdgeFaceStencil[0];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil1 = fd_vertexToDiagonalEdgeFaceStencil[1];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil2 = fd_vertexToDiagonalEdgeFaceStencil[2];
  const double asdf_fd_vertexToDiagonalEdgeFaceStencil3 = fd_vertexToDiagonalEdgeFaceStencil[3];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil0 = fd_vertexToVerticalEdgeFaceStencil[0];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil1 = fd_vertexToVerticalEdgeFaceStencil[1];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil2 = fd_vertexToVerticalEdgeFaceStencil[2];
  const double asdf_fd_vertexToVerticalEdgeFaceStencil3 = fd_vertexToVerticalEdgeFaceStencil[3];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil0 = fd_vertexToHorizontalEdgeFaceStencil[0];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil1 = fd_vertexToHorizontalEdgeFaceStencil[1];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil2 = fd_vertexToHorizontalEdgeFaceStencil[2];
  const double asdf_fd_vertexToHorizontalEdgeFaceStencil3 = fd_vertexToHorizontalEdgeFaceStencil[3];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToDiagonalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToDiagonalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16387] + asdf_fd_vertexToDiagonalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToDiagonalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToVerticalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToVerticalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToVerticalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + asdf_fd_vertexToVerticalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
  for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = asdf_fd_vertexToHorizontalEdgeFaceStencil0*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + asdf_fd_vertexToHorizontalEdgeFaceStencil1*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + asdf_fd_vertexToHorizontalEdgeFaceStencil2*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + asdf_fd_vertexToHorizontalEdgeFaceStencil3*fd_vertexFaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16386] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)];
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

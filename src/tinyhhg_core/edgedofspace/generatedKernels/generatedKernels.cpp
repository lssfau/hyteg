
#include "generatedKernels.hpp"
#include "core/logging/Logging.h"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {


static void apply_2D_macroface_edgedof_to_edgedof_replace_level_2(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_3(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_4(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_5(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_6(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_7(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_8(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_9(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_10(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_11(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_12(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_13(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_14(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2) - 1];
    }
  }
  for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 - 1)) / 2) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / 2) - ((ctr_2*(ctr_2 + 1)) / 2)];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_any(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil, int64_t level)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2) - 1];
    }
    for (int ctr_1 = (1 << level) - 1; ctr_1 < (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2) - 1];
    }
  }
  for (int ctr_2 = 1; ctr_2 < (1 << level) - 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2)] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2) - 1];
    }
    for (int ctr_1 = -ctr_2 + (1 << level) - 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / 2) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2) - 1];
    }
  }
  for (int ctr_2 = (1 << level) - 1; ctr_2 < (1 << level); ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2)] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + ((((1 << level) + 1)*(1 << level)) / 2)] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / 2) + 2*((((1 << level) + 1)*(1 << level)) / 2)];
    }
  }
}





static void apply_2D_macroface_edgedof_to_edgedof_replace(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil, int64_t level)
{
  switch( level )
  {
    case 2:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_2(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 3:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_3(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 4:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_4(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 5:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_5(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 6:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_6(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 7:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_7(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 8:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_8(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 9:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_9(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 10:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_10(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 11:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_11(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 12:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_12(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 13:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_13(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 14:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_14(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    default:
      apply_2D_macroface_edgedof_to_edgedof_replace_level_any(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil, level);
      break;
  }
}


static void apply_2D_macroface_edgedof_to_edgedof_add_level_2(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_3(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_4(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_5(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_6(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_7(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_8(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_9(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_10(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_11(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_12(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_13(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_14(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
  for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_add_level_any(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil, int64_t level)
{
  const double tmpconst_fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double tmpconst_fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double tmpconst_fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double tmpconst_fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double tmpconst_fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double tmpconst_fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double tmpconst_fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double tmpconst_fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double tmpconst_fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double tmpconst_fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  const double tmpconst_fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double tmpconst_fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double tmpconst_fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double tmpconst_fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double tmpconst_fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
    }
    for (int ctr_1 = (1 << level) - 1; ctr_1 < (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
    }
  }
  for (int ctr_2 = 1; ctr_2 < (1 << level) - 1; ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
    }
    for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level) - 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + tmpconst_fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
    }
    for (int ctr_1 = -ctr_2 + (1 << level) - 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = tmpconst_fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + tmpconst_fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
    }
  }
  for (int ctr_2 = (1 << level) - 1; ctr_2 < (1 << level); ctr_2 += 1)
  {
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = tmpconst_fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + tmpconst_fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + tmpconst_fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + tmpconst_fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
    }
  }
}




static void apply_2D_macroface_edgedof_to_edgedof_add(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil, int64_t level)
{
  switch( level )
  {
    case 2:
      apply_2D_macroface_edgedof_to_edgedof_add_level_2(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 3:
      apply_2D_macroface_edgedof_to_edgedof_add_level_3(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 4:
      apply_2D_macroface_edgedof_to_edgedof_add_level_4(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 5:
      apply_2D_macroface_edgedof_to_edgedof_add_level_5(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 6:
      apply_2D_macroface_edgedof_to_edgedof_add_level_6(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 7:
      apply_2D_macroface_edgedof_to_edgedof_add_level_7(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 8:
      apply_2D_macroface_edgedof_to_edgedof_add_level_8(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 9:
      apply_2D_macroface_edgedof_to_edgedof_add_level_9(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 10:
      apply_2D_macroface_edgedof_to_edgedof_add_level_10(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 11:
      apply_2D_macroface_edgedof_to_edgedof_add_level_11(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 12:
      apply_2D_macroface_edgedof_to_edgedof_add_level_12(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 13:
      apply_2D_macroface_edgedof_to_edgedof_add_level_13(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    case 14:
      apply_2D_macroface_edgedof_to_edgedof_add_level_14(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil);
      break;
    default:
      apply_2D_macroface_edgedof_to_edgedof_add_level_any(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil, level);
      break;
  }
}





void applyReplace( double*                fd_edgeFaceDst,
                   double*          fd_edgeFaceSrc,
                   double*          fd_edgeFaceStencil,
                   walberla::uint_t level )
{
  apply_2D_macroface_edgedof_to_edgedof_replace(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil, static_cast< int64_t >(level));
}

void applyAdd( double* fd_edgeFaceDst, double* fd_edgeFaceSrc, double* fd_edgeFaceStencil, walberla::uint_t level )
{
  apply_2D_macroface_edgedof_to_edgedof_add(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil, static_cast< int64_t >(level));

}


} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
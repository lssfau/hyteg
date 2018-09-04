
#include "generatedKernels.hpp"
#include "core/logging/Logging.h"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_edgedof_to_edgedof_replace_level_2(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (20 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (20 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 5] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2)];
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (20 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(20 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (20 / 2) + 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 4] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 8] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 3];
        fd_edgeFaceDst[-(0 / 2) + 2*(20 / 2) + 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 2];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 5] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 5] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (20 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 5];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 5] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 5] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 5];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 6] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 6] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 6] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 6] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6];
          fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 6] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 6] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 7] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 7] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 7] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 7] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 8] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 12] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7];
          fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6];
        }
      }
      {
        {
          fd_edgeFaceDst[-(6 / 2) + 10] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 10];
          fd_edgeFaceDst[-(6 / 2) + (20 / 2) + 10] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 15] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 10];
        }
        {
          fd_edgeFaceDst[-(6 / 2) + 11] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11];
          fd_edgeFaceDst[-(6 / 2) + (20 / 2) + 11] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 12] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 16] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11];
          fd_edgeFaceDst[-(6 / 2) + 2*(20 / 2) + 11] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (12 / 2) + 15] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 15] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (20 / 2) + 15] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(20 / 2) + 15];
      fd_edgeFaceDst[ctr_1 - (12 / 2) + (20 / 2) + 15] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (20 / 2) + 15] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 15] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(20 / 2) + 16] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 20] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(20 / 2) + 15];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_3(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (72 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (72 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 9] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2)];
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (72 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (72 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(72 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 9] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(72 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(72 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(72 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (72 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (72 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (72 / 2) + 7] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 7] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 7] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 8] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 16] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 7];
        fd_edgeFaceDst[-(0 / 2) + 2*(72 / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 6];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeFaceStencil2*fd_edgeFaceSrc[9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8] + fd_edgeFaceStencil3*fd_edgeFaceSrc[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_edgeFaceStencil9*fd_edgeFaceSrc[9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_edgeFaceStencil0*fd_edgeFaceSrc[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil1*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil4*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7];
        fd_edgeFaceDst[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_edgeFaceStencil5*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil6*fd_edgeFaceSrc[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil7*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8] + fd_edgeFaceStencil8*fd_edgeFaceSrc[8*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16] + fd_edgeFaceStencil9*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7];
        fd_edgeFaceDst[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[8*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 6];
      }
    }
    {
      {
        fd_edgeFaceDst[-(56 / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 63];
        fd_edgeFaceDst[-(56 / 2) + (72 / 2) + 63] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 64] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(72 / 2) + 72] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 63];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (56 / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (56 / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (56 / 2) + 2*(72 / 2) + 63];
      }
      {
        fd_edgeFaceDst[-(56 / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 63];
        fd_edgeFaceDst[-(56 / 2) + (72 / 2) + 63] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 64] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(72 / 2) + 72] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(56 / 2) + 2*(72 / 2) + 63];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_4(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (272 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (272 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 17] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2)];
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (272 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(272 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (272 / 2) + 15] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 15] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 15] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 16] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 32] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 15];
        fd_edgeFaceDst[-(0 / 2) + 2*(272 / 2) + 15] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 15] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 15] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 15] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 31] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 14];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeFaceStencil2*fd_edgeFaceSrc[17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16] + fd_edgeFaceStencil3*fd_edgeFaceSrc[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_edgeFaceStencil9*fd_edgeFaceSrc[17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_edgeFaceStencil0*fd_edgeFaceSrc[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil1*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil4*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15];
        fd_edgeFaceDst[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_edgeFaceStencil5*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil6*fd_edgeFaceSrc[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil7*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16] + fd_edgeFaceStencil8*fd_edgeFaceSrc[16*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32] + fd_edgeFaceStencil9*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15];
        fd_edgeFaceDst[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_edgeFaceStencil10*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil11*fd_edgeFaceSrc[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil12*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil13*fd_edgeFaceSrc[16*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 31] + fd_edgeFaceStencil14*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 14];
      }
    }
    {
      {
        fd_edgeFaceDst[-(240 / 2) + 255] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 255];
        fd_edgeFaceDst[-(240 / 2) + (272 / 2) + 255] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 256] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(272 / 2) + 272] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 255];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (240 / 2) + 255] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (240 / 2) + 255] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (240 / 2) + 2*(272 / 2) + 255];
      }
      {
        fd_edgeFaceDst[-(240 / 2) + 255] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 255];
        fd_edgeFaceDst[-(240 / 2) + (272 / 2) + 255] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 256] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(272 / 2) + 272] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(240 / 2) + 2*(272 / 2) + 255];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_5(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1056 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 33] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2)];
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1056 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1056 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 32] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (1056 / 2) + 31] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 31] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 31] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 32] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 64] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 31];
        fd_edgeFaceDst[-(0 / 2) + 2*(1056 / 2) + 31] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 31] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 31] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 31] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 63] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 30];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 33] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 33] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 33] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 33];
          fd_edgeFaceDst[-(2 / 2) + (1056 / 2) + 33] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 33] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 33] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 34] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 66] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 33];
        }
        for (int ctr_1 = 1; ctr_1 < 30; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 33] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (1056 / 2) + 33] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 34] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 66] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 65] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 32];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 30] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 31] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 63];
          fd_edgeFaceDst[-(2 / 2) + (1056 / 2) + 63] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 63] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 63] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 64] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 96] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 63];
          fd_edgeFaceDst[-(2 / 2) + 2*(1056 / 2) + 63] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(1056 / 2) + 63] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 63] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 63] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 95] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (1056 / 2) + 62];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 30; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeFaceStencil2*fd_edgeFaceSrc[33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 32] + fd_edgeFaceStencil3*fd_edgeFaceSrc[33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_edgeFaceStencil9*fd_edgeFaceSrc[33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 32] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
        }
        {
          fd_edgeFaceDst[32*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] = fd_edgeFaceStencil0*fd_edgeFaceSrc[32*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil1*fd_edgeFaceSrc[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil4*fd_edgeFaceSrc[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31];
          fd_edgeFaceDst[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] = fd_edgeFaceStencil5*fd_edgeFaceSrc[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil6*fd_edgeFaceSrc[32*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil7*fd_edgeFaceSrc[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 32] + fd_edgeFaceStencil8*fd_edgeFaceSrc[32*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 64] + fd_edgeFaceStencil9*fd_edgeFaceSrc[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31];
          fd_edgeFaceDst[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] = fd_edgeFaceStencil10*fd_edgeFaceSrc[32*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil11*fd_edgeFaceSrc[32*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil12*fd_edgeFaceSrc[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 31] + fd_edgeFaceStencil13*fd_edgeFaceSrc[32*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 63] + fd_edgeFaceStencil14*fd_edgeFaceSrc[32*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 30];
        }
      }
      {
        {
          fd_edgeFaceDst[-(930 / 2) + 990] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(930 / 2) + 990] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(870 / 2) + (1056 / 2) + 957] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(870 / 2) + 2*(1056 / 2) + 958] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 990];
          fd_edgeFaceDst[-(930 / 2) + (1056 / 2) + 990] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(930 / 2) + 990] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(992 / 2) + 1023] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 990];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (930 / 2) + 990] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 990] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (870 / 2) + (1056 / 2) + 957] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 2*(1056 / 2) + 958] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 990];
          fd_edgeFaceDst[ctr_1 - (930 / 2) + (1056 / 2) + 990] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 990] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1023] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 990];
        }
        {
          fd_edgeFaceDst[-(930 / 2) + 991] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(930 / 2) + 991] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(870 / 2) + (1056 / 2) + 958] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(870 / 2) + 2*(1056 / 2) + 959] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 991] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991];
          fd_edgeFaceDst[-(930 / 2) + (1056 / 2) + 991] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 991] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(930 / 2) + 991] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 992] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(992 / 2) + 1024] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991];
          fd_edgeFaceDst[-(930 / 2) + 2*(1056 / 2) + 991] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(930 / 2) + 991] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 991] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(992 / 2) + 1023] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 990];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (992 / 2) + 1023] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1023] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (992 / 2) + (1056 / 2) + 1023] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 2*(1056 / 2) + 1023];
      fd_edgeFaceDst[ctr_1 - (992 / 2) + (1056 / 2) + 1023] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (992 / 2) + (1056 / 2) + 1023] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1023] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 2*(1056 / 2) + 1024] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (1056 / 2) + 1056] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 2*(1056 / 2) + 1023];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_6(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4160 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 65] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2)];
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4160 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4160 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4160 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 65] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4160 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4160 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4160 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4160 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 64] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4160 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (4160 / 2) + 63] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2) + 63] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 63] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2) + 64] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 128] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2) + 63];
        fd_edgeFaceDst[-(0 / 2) + 2*(4160 / 2) + 63] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2) + 63] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 63] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2) + 63] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 127] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2) + 62];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeFaceStencil2*fd_edgeFaceSrc[65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 64] + fd_edgeFaceStencil3*fd_edgeFaceSrc[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_edgeFaceStencil9*fd_edgeFaceSrc[65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 64] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 64] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63];
        fd_edgeFaceDst[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_edgeFaceStencil5*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil6*fd_edgeFaceSrc[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil7*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 64] + fd_edgeFaceStencil8*fd_edgeFaceSrc[64*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128] + fd_edgeFaceStencil9*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63];
        fd_edgeFaceDst[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_edgeFaceStencil10*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil11*fd_edgeFaceSrc[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil12*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil13*fd_edgeFaceSrc[64*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 127] + fd_edgeFaceStencil14*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 62];
      }
    }
    {
      {
        fd_edgeFaceDst[-(4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
        fd_edgeFaceDst[-(4032 / 2) + (4160 / 2) + 4095] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4096] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(4160 / 2) + 4160] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + 2*(4160 / 2) + 4095];
      }
      {
        fd_edgeFaceDst[-(4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
        fd_edgeFaceDst[-(4032 / 2) + (4160 / 2) + 4095] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4096] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(4160 / 2) + 4160] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_7(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16512 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 129] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2)];
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16512 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16512 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16512 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 129] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16512 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16512 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16512 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16512 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 128] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16512 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (16512 / 2) + 127] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2) + 127] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 127] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2) + 128] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 256] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2) + 127];
        fd_edgeFaceDst[-(0 / 2) + 2*(16512 / 2) + 127] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2) + 127] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 127] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2) + 127] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 255] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2) + 126];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeFaceStencil2*fd_edgeFaceSrc[129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 128] + fd_edgeFaceStencil3*fd_edgeFaceSrc[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_edgeFaceStencil9*fd_edgeFaceSrc[129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 128] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_edgeFaceStencil0*fd_edgeFaceSrc[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil1*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil4*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127];
        fd_edgeFaceDst[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_edgeFaceStencil5*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil6*fd_edgeFaceSrc[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil7*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 128] + fd_edgeFaceStencil8*fd_edgeFaceSrc[128*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256] + fd_edgeFaceStencil9*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127];
        fd_edgeFaceDst[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_edgeFaceStencil10*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil11*fd_edgeFaceSrc[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil12*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil13*fd_edgeFaceSrc[128*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 255] + fd_edgeFaceStencil14*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 126];
      }
    }
    {
      {
        fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
        fd_edgeFaceDst[-(16256 / 2) + (16512 / 2) + 16383] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16384] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(16512 / 2) + 16512] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + 2*(16512 / 2) + 16383];
      }
      {
        fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
        fd_edgeFaceDst[-(16256 / 2) + (16512 / 2) + 16383] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16384] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(16512 / 2) + 16512] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_8(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (65792 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 257] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2)];
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (65792 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (65792 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(65792 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 257] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(65792 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(65792 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(65792 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (65792 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 256] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (65792 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (65792 / 2) + 255] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 255] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 255] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 256] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 512] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 255];
        fd_edgeFaceDst[-(0 / 2) + 2*(65792 / 2) + 255] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 255] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 255] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 255] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 254];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 257] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 257] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 257];
          fd_edgeFaceDst[-(2 / 2) + (65792 / 2) + 257] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 257] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 258] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 514] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 257];
        }
        for (int ctr_1 = 1; ctr_1 < 254; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 257] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 257] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (65792 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(65792 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(65792 / 2) + 257];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (65792 / 2) + 257] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 257] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(65792 / 2) + 258] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 514] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(65792 / 2) + 257];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(65792 / 2) + 257] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(65792 / 2) + 257] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 257] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 513] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (65792 / 2) + 256];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 511] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 511] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 511];
          fd_edgeFaceDst[-(2 / 2) + (65792 / 2) + 511] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 511] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 512] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 768] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 511];
          fd_edgeFaceDst[-(2 / 2) + 2*(65792 / 2) + 511] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 511] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 511] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 510];
        }
      }
      {
        for (int ctr_1 = 0; ctr_1 < 254; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (6 / 2) + 514] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 514] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (65792 / 2) + 257] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(65792 / 2) + 258] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (65792 / 2) + 514] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(65792 / 2) + 514];
          fd_edgeFaceDst[ctr_1 - (6 / 2) + (65792 / 2) + 514] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (65792 / 2) + 514] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 514] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(65792 / 2) + 515] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 771] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(65792 / 2) + 514];
          if (ctr_1 > 0)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(65792 / 2) + 514] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(65792 / 2) + 514] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 514] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (65792 / 2) + 514] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 770] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (65792 / 2) + 513];
          }
        }
        for (int ctr_2 = 3; ctr_2 < 253; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 256] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            if (ctr_1 > 0)
            {
              fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (64262 / 2) + 65021] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 65021] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (63756 / 2) + (65792 / 2) + 64764] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (63756 / 2) + 2*(65792 / 2) + 64765] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65021];
          fd_edgeFaceDst[ctr_1 - (64262 / 2) + (65792 / 2) + 65021] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 65021] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65022] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 65278] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65021];
          if (ctr_1 > 0)
          {
            fd_edgeFaceDst[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65021] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65021] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 65021] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 65277] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + (65792 / 2) + 65020];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(64770 / 2) + 65278] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(64770 / 2) + 65278] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(64262 / 2) + 2*(65792 / 2) + 65022] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65278];
          fd_edgeFaceDst[-(64770 / 2) + (65792 / 2) + 65278] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(64770 / 2) + 65278] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(65280 / 2) + 65535] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65278];
        }
        {
          fd_edgeFaceDst[-(64770 / 2) + 65279] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(64770 / 2) + 65279] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(64262 / 2) + (65792 / 2) + 65022] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(64262 / 2) + 2*(65792 / 2) + 65023] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65279] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279];
          fd_edgeFaceDst[-(64770 / 2) + (65792 / 2) + 65279] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65279] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(64770 / 2) + 65279] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65280] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(65280 / 2) + 65536] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279];
          fd_edgeFaceDst[-(64770 / 2) + 2*(65792 / 2) + 65279] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(64770 / 2) + 65279] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65279] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(65280 / 2) + 65535] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (65280 / 2) + 65535] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 65535] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + (65792 / 2) + 65535] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 2*(65792 / 2) + 65535];
      fd_edgeFaceDst[ctr_1 - (65280 / 2) + (65792 / 2) + 65535] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + (65792 / 2) + 65535] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 65535] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 2*(65792 / 2) + 65536] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (65792 / 2) + 65792] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 2*(65792 / 2) + 65535];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_9(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (262656 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 513] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2)];
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (262656 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(262656 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 512] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (262656 / 2) + 511] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 511] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 511] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 512] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 1024] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 511];
        fd_edgeFaceDst[-(0 / 2) + 2*(262656 / 2) + 511] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 511] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 511] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 511] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 1023] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 510];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeFaceStencil2*fd_edgeFaceSrc[513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 512] + fd_edgeFaceStencil3*fd_edgeFaceSrc[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_edgeFaceStencil9*fd_edgeFaceSrc[513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 512] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_edgeFaceStencil0*fd_edgeFaceSrc[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil1*fd_edgeFaceSrc[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil4*fd_edgeFaceSrc[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511];
        fd_edgeFaceDst[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_edgeFaceStencil5*fd_edgeFaceSrc[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil6*fd_edgeFaceSrc[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil7*fd_edgeFaceSrc[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 512] + fd_edgeFaceStencil8*fd_edgeFaceSrc[512*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024] + fd_edgeFaceStencil9*fd_edgeFaceSrc[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511];
        fd_edgeFaceDst[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] = fd_edgeFaceStencil10*fd_edgeFaceSrc[512*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil11*fd_edgeFaceSrc[512*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil12*fd_edgeFaceSrc[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 511] + fd_edgeFaceStencil13*fd_edgeFaceSrc[512*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1023] + fd_edgeFaceStencil14*fd_edgeFaceSrc[512*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 510];
      }
    }
    {
      {
        fd_edgeFaceDst[-(261632 / 2) + 262143] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262143];
        fd_edgeFaceDst[-(261632 / 2) + (262656 / 2) + 262143] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262144] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(262656 / 2) + 262656] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262143];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (261632 / 2) + 262143] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + 262143] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + 2*(262656 / 2) + 262143];
      }
      {
        fd_edgeFaceDst[-(261632 / 2) + 262143] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262143];
        fd_edgeFaceDst[-(261632 / 2) + (262656 / 2) + 262143] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262144] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(262656 / 2) + 262656] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(261632 / 2) + 2*(262656 / 2) + 262143];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_10(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (1049600 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 1025] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2)];
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (1049600 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(1049600 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1024] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (1049600 / 2) + 1023] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1023] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 1023] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1024] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2048] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1023];
        fd_edgeFaceDst[-(0 / 2) + 2*(1049600 / 2) + 1023] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1023] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 1023] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1023] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1022];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 1025] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 1025] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 1025];
          fd_edgeFaceDst[-(2 / 2) + (1049600 / 2) + 1025] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 1025] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 1026] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2050] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 1025];
        }
        for (int ctr_1 = 1; ctr_1 < 1022; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 1025] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1026] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2049] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1024];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 2047] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1022] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1023] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2047] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2047];
          fd_edgeFaceDst[-(2 / 2) + (1049600 / 2) + 2047] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2047] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2048] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 3072] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2047];
          fd_edgeFaceDst[-(2 / 2) + 2*(1049600 / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 3071] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2046];
        }
      }
      {
        for (int ctr_1 = 0; ctr_1 < 1022; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (6 / 2) + 2050] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1026] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050];
          fd_edgeFaceDst[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2051] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3075] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050];
          if (ctr_1 > 0)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3074] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2049];
          }
        }
        for (int ctr_2 = 3; ctr_2 < 1021; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1024] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            if (ctr_1 > 0)
            {
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (1043462 / 2) + 1046525] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 1046525] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045501] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 2*(1049600 / 2) + 1046525];
          fd_edgeFaceDst[ctr_1 - (1043462 / 2) + (1049600 / 2) + 1046525] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 1046525] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 2*(1049600 / 2) + 1046526] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + 1047550] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 2*(1049600 / 2) + 1046525];
          if (ctr_1 > 0)
          {
            fd_edgeFaceDst[ctr_1 - (1043462 / 2) + 2*(1049600 / 2) + 1046525] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 2*(1049600 / 2) + 1046525] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 1046525] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + 1047549] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + (1049600 / 2) + 1046524];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(1045506 / 2) + 1047550] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1045506 / 2) + 1047550] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046526] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047550] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047550];
          fd_edgeFaceDst[-(1045506 / 2) + (1049600 / 2) + 1047550] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047550] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(1045506 / 2) + 1047550] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(1047552 / 2) + 1048575] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047550];
        }
        {
          fd_edgeFaceDst[-(1045506 / 2) + 1047551] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046526] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046527] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047551] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047551];
          fd_edgeFaceDst[-(1045506 / 2) + (1049600 / 2) + 1047551] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047551] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047552] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(1047552 / 2) + 1048576] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047551];
          fd_edgeFaceDst[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047551] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(1047552 / 2) + 1048575] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047550];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (1047552 / 2) + 1048575] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 1048575] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + (1049600 / 2) + 1047550] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + 2*(1049600 / 2) + 1047551] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + (1049600 / 2) + 1048575] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 2*(1049600 / 2) + 1048575];
      fd_edgeFaceDst[ctr_1 - (1047552 / 2) + (1049600 / 2) + 1048575] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + (1049600 / 2) + 1048575] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 1048575] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 2*(1049600 / 2) + 1048576] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (1049600 / 2) + 1049600] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 2*(1049600 / 2) + 1048575];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_11(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (4196352 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2049] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2)];
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (4196352 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4196352 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4196352 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2049] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4196352 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(4196352 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4196352 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4196352 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2048] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4196352 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (4196352 / 2) + 2047] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2) + 2047] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 2047] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2) + 2048] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 4096] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2) + 2047];
        fd_edgeFaceDst[-(0 / 2) + 2*(4196352 / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 4095] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2) + 2046];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeFaceStencil2*fd_edgeFaceSrc[2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2048] + fd_edgeFaceStencil3*fd_edgeFaceSrc[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_edgeFaceStencil9*fd_edgeFaceSrc[2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2048] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2048] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_edgeFaceStencil0*fd_edgeFaceSrc[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil1*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil4*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047];
        fd_edgeFaceDst[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_edgeFaceStencil5*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil6*fd_edgeFaceSrc[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil7*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2048] + fd_edgeFaceStencil8*fd_edgeFaceSrc[2048*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096] + fd_edgeFaceStencil9*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047];
        fd_edgeFaceDst[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[2048*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4095] + fd_edgeFaceStencil14*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2046];
      }
    }
    {
      {
        fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
        fd_edgeFaceDst[-(4192256 / 2) + (4196352 / 2) + 4194303] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194304] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(4196352 / 2) + 4196352] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + 2*(4196352 / 2) + 4194303];
      }
      {
        fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
        fd_edgeFaceDst[-(4192256 / 2) + (4196352 / 2) + 4194303] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194304] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(4196352 / 2) + 4196352] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_12(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (16781312 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 4097] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2)];
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (16781312 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(16781312 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4096] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (16781312 / 2) + 4095] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4095] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 4095] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 4096] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 8192] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 4095];
        fd_edgeFaceDst[-(0 / 2) + 2*(16781312 / 2) + 4095] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 4095] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 4095] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4095] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 8191] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4094];
      }
    }
    for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
    {
      {
        fd_edgeFaceDst[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeFaceStencil2*fd_edgeFaceSrc[4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4096] + fd_edgeFaceStencil3*fd_edgeFaceSrc[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_edgeFaceStencil9*fd_edgeFaceSrc[4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4096] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        fd_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
      }
      {
        fd_edgeFaceDst[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095];
        fd_edgeFaceDst[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_edgeFaceStencil5*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil6*fd_edgeFaceSrc[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil7*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4096] + fd_edgeFaceStencil8*fd_edgeFaceSrc[4096*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192] + fd_edgeFaceStencil9*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095];
        fd_edgeFaceDst[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_edgeFaceStencil10*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil11*fd_edgeFaceSrc[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil12*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil13*fd_edgeFaceSrc[4096*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8191] + fd_edgeFaceStencil14*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4094];
      }
    }
    {
      {
        fd_edgeFaceDst[-(16773120 / 2) + 16777215] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777215];
        fd_edgeFaceDst[-(16773120 / 2) + (16781312 / 2) + 16777215] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777216] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(16781312 / 2) + 16781312] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777215];
      }
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (16773120 / 2) + 16777215] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + 16777215] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + 2*(16781312 / 2) + 16777215];
      }
      {
        fd_edgeFaceDst[-(16773120 / 2) + 16777215] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777215];
        fd_edgeFaceDst[-(16773120 / 2) + (16781312 / 2) + 16777215] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777216] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(16781312 / 2) + 16781312] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(16773120 / 2) + 2*(16781312 / 2) + 16777215];
      }
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_13(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (67117056 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 8193] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2)];
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (67117056 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(67117056 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8192] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (67117056 / 2) + 8191] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8191] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 8191] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 8192] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 16384] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 8191];
        fd_edgeFaceDst[-(0 / 2) + 2*(67117056 / 2) + 8191] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 8191] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 8191] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8191] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 16383] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8190];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 8193] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 8193] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 8193];
          fd_edgeFaceDst[-(2 / 2) + (67117056 / 2) + 8193] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 8193] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 8194] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 16386] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 8193];
        }
        for (int ctr_1 = 1; ctr_1 < 8190; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 8193] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8194] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16386] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16385] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8192];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8190] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 8191] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 16383];
          fd_edgeFaceDst[-(2 / 2) + (67117056 / 2) + 16383] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 16383] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 16383] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 16384] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 24576] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 16383];
          fd_edgeFaceDst[-(2 / 2) + 2*(67117056 / 2) + 16383] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 16383] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 16383] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 16383] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 24575] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 16382];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 16386] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 16386] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 8194] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 16386] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 16386];
            fd_edgeFaceDst[-(6 / 2) + (67117056 / 2) + 16386] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 16386] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 16386] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 16387] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 24579] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 16386];
          }
          for (int ctr_1 = 1; ctr_1 < 8189; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 16386] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16386] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8194] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16386] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16387] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 24579] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16386] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 24578] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16385];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 24575] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 24575] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (67117056 / 2) + 16382] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(67117056 / 2) + 16383] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 24575] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 24575];
            fd_edgeFaceDst[-(6 / 2) + (67117056 / 2) + 24575] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 24575] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 24575] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 24576] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 32768] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 24575];
            fd_edgeFaceDst[-(6 / 2) + 2*(67117056 / 2) + 24575] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(67117056 / 2) + 24575] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 24575] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 24575] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 32767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (67117056 / 2) + 24574];
          }
        }
        for (int ctr_2 = 3; ctr_2 < 8189; ctr_2 += 1)
        {
          {
            fd_edgeFaceDst[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeFaceStencil2*fd_edgeFaceSrc[8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8192] + fd_edgeFaceStencil3*fd_edgeFaceSrc[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_edgeFaceStencil9*fd_edgeFaceSrc[8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          }
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8192] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            fd_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
          }
          {
            fd_edgeFaceDst[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_edgeFaceStencil0*fd_edgeFaceSrc[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil1*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil4*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191];
            fd_edgeFaceDst[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_edgeFaceStencil5*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil6*fd_edgeFaceSrc[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil7*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8192] + fd_edgeFaceStencil8*fd_edgeFaceSrc[8192*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384] + fd_edgeFaceStencil9*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191];
            fd_edgeFaceDst[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_edgeFaceStencil10*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil11*fd_edgeFaceSrc[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil12*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil13*fd_edgeFaceSrc[8192*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16383] + fd_edgeFaceStencil14*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8190];
          }
        }
        {
          {
            fd_edgeFaceDst[-(67067910 / 2) + 67092477] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67067910 / 2) + 67092477] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67051532 / 2) + (67117056 / 2) + 67084284] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67051532 / 2) + 2*(67117056 / 2) + 67084285] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092477];
            fd_edgeFaceDst[-(67067910 / 2) + (67117056 / 2) + 67092477] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(67067910 / 2) + 67092477] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092478] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(67084290 / 2) + 67100670] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092477];
          }
          for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (67067910 / 2) + 67092477] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 67092477] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (67051532 / 2) + (67117056 / 2) + 67084284] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (67051532 / 2) + 2*(67117056 / 2) + 67084285] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092477];
            fd_edgeFaceDst[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092477] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 67092477] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092478] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 67100670] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092477];
            fd_edgeFaceDst[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092477] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092477] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 67092477] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 67100669] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092476];
          }
          {
            fd_edgeFaceDst[-(67067910 / 2) + 67092479] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67067910 / 2) + 67092479] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67051532 / 2) + (67117056 / 2) + 67084286] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67051532 / 2) + 2*(67117056 / 2) + 67084287] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092479] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092479];
            fd_edgeFaceDst[-(67067910 / 2) + (67117056 / 2) + 67092479] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092479] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(67067910 / 2) + 67092479] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092480] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(67084290 / 2) + 67100672] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092479];
            fd_edgeFaceDst[-(67067910 / 2) + 2*(67117056 / 2) + 67092479] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092479] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(67067910 / 2) + 67092479] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092479] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(67084290 / 2) + 67100671] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092478];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(67084290 / 2) + 67100670] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67084290 / 2) + 67100670] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092478] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100670];
          fd_edgeFaceDst[-(67084290 / 2) + (67117056 / 2) + 67100670] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(67084290 / 2) + 67100670] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(67100672 / 2) + 67108863] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100670];
        }
        {
          fd_edgeFaceDst[-(67084290 / 2) + 67100671] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67084290 / 2) + 67100671] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67067910 / 2) + (67117056 / 2) + 67092478] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67067910 / 2) + 2*(67117056 / 2) + 67092479] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100671] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671];
          fd_edgeFaceDst[-(67084290 / 2) + (67117056 / 2) + 67100671] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100671] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(67084290 / 2) + 67100671] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100672] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(67100672 / 2) + 67108864] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671];
          fd_edgeFaceDst[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(67084290 / 2) + 67100671] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100671] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(67100672 / 2) + 67108863] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100670];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (67100672 / 2) + 67108863] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 67108863] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + (67117056 / 2) + 67108863] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 2*(67117056 / 2) + 67108863];
      fd_edgeFaceDst[ctr_1 - (67100672 / 2) + (67117056 / 2) + 67108863] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + (67117056 / 2) + 67108863] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 67108863] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 2*(67117056 / 2) + 67108864] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (67117056 / 2) + 67117056] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 2*(67117056 / 2) + 67108863];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_14(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      fd_edgeFaceDst[-(0 / 2) + (268451840 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 16385] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2)];
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (0 / 2) + (268451840 / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2)];
        fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(268451840 / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16384] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2) - 1];
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (268451840 / 2) + 16383] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16383] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 16383] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 16384] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 32768] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 16383];
        fd_edgeFaceDst[-(0 / 2) + 2*(268451840 / 2) + 16383] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 16383] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 16383] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16383] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 32767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16382];
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 16385] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 16385] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 16385] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 16385];
          fd_edgeFaceDst[-(2 / 2) + (268451840 / 2) + 16385] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 16385] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 16385] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 16386] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 32770] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 16385];
        }
        for (int ctr_1 = 1; ctr_1 < 16382; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 16385] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16386] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 32770] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 32769] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16384];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 32767] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 32767] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16382] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 16383] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 32767] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 32767];
          fd_edgeFaceDst[-(2 / 2) + (268451840 / 2) + 32767] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 32767] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 32767] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 32768] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 49152] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 32767];
          fd_edgeFaceDst[-(2 / 2) + 2*(268451840 / 2) + 32767] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(268451840 / 2) + 32767] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 32767] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 32767] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 49151] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (268451840 / 2) + 32766];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 16382; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeFaceStencil2*fd_edgeFaceSrc[16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16384] + fd_edgeFaceStencil3*fd_edgeFaceSrc[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_edgeFaceStencil9*fd_edgeFaceSrc[16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16384] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
        }
        {
          fd_edgeFaceDst[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383];
          fd_edgeFaceDst[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_edgeFaceStencil5*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil6*fd_edgeFaceSrc[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil7*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16384] + fd_edgeFaceStencil8*fd_edgeFaceSrc[16384*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32768] + fd_edgeFaceStencil9*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383];
          fd_edgeFaceDst[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_edgeFaceStencil10*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil11*fd_edgeFaceSrc[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil12*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil13*fd_edgeFaceSrc[16384*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16382];
        }
      }
      {
        {
          fd_edgeFaceDst[-(268386306 / 2) + 268419070] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(268386306 / 2) + 268419070] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(268353542 / 2) + (268451840 / 2) + 268402685] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(268353542 / 2) + 2*(268451840 / 2) + 268402686] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419070];
          fd_edgeFaceDst[-(268386306 / 2) + (268451840 / 2) + 268419070] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(268386306 / 2) + 268419070] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(268419072 / 2) + 268435455] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419070];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (268386306 / 2) + 268419070] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 268419070] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + (268451840 / 2) + 268402685] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + 2*(268451840 / 2) + 268402686] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419070];
          fd_edgeFaceDst[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 268419070] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 268435455] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419070];
        }
        {
          fd_edgeFaceDst[-(268386306 / 2) + 268419071] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(268386306 / 2) + 268419071] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(268353542 / 2) + (268451840 / 2) + 268402686] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(268353542 / 2) + 2*(268451840 / 2) + 268402687] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419071] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071];
          fd_edgeFaceDst[-(268386306 / 2) + (268451840 / 2) + 268419071] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419071] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(268386306 / 2) + 268419071] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419072] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(268419072 / 2) + 268435456] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071];
          fd_edgeFaceDst[-(268386306 / 2) + 2*(268451840 / 2) + 268419071] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(268386306 / 2) + 268419071] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419071] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(268419072 / 2) + 268435455] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419070];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (268419072 / 2) + 268435455] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 268435455] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + (268451840 / 2) + 268435455] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 2*(268451840 / 2) + 268435455];
      fd_edgeFaceDst[ctr_1 - (268419072 / 2) + (268451840 / 2) + 268435455] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + (268451840 / 2) + 268435455] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 268435455] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 2*(268451840 / 2) + 268435456] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (268451840 / 2) + 268451840] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 2*(268451840 / 2) + 268435455];
    }
  }
}



static void apply_2D_macroface_edgedof_to_edgedof_replace_level_any(double * fd_edgeFaceDst, double * fd_edgeFaceSrc, double * fd_edgeFaceStencil, int64_t level)
{
  const double fd_edgeFaceStencil0 = fd_edgeFaceStencil[0];
  const double fd_edgeFaceStencil1 = fd_edgeFaceStencil[1];
  const double fd_edgeFaceStencil2 = fd_edgeFaceStencil[2];
  const double fd_edgeFaceStencil3 = fd_edgeFaceStencil[3];
  const double fd_edgeFaceStencil4 = fd_edgeFaceStencil[4];
  const double fd_edgeFaceStencil5 = fd_edgeFaceStencil[5];
  const double fd_edgeFaceStencil6 = fd_edgeFaceStencil[6];
  const double fd_edgeFaceStencil7 = fd_edgeFaceStencil[7];
  const double fd_edgeFaceStencil8 = fd_edgeFaceStencil[8];
  const double fd_edgeFaceStencil9 = fd_edgeFaceStencil[9];
  const double fd_edgeFaceStencil10 = fd_edgeFaceStencil[10];
  const double fd_edgeFaceStencil11 = fd_edgeFaceStencil[11];
  const double fd_edgeFaceStencil12 = fd_edgeFaceStencil[12];
  const double fd_edgeFaceStencil13 = fd_edgeFaceStencil[13];
  const double fd_edgeFaceStencil14 = fd_edgeFaceStencil[14];
  {
    {
      if ((1 << level) > 0)
      {
        fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
      {
        {
          if ((1 << level) > 1)
          {
            fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
          }
          fd_edgeFaceDst[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2)];
        }
        {
          {
            if ((1 << level) > 2)
            {
              fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 3] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
            }
            fd_edgeFaceDst[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
          }
          for (int ctr_1 = 3; ctr_1 < (1 << level) - 3; ctr_1 += 1)
          {
            if (ctr_1 < (1 << level))
            {
              fd_edgeFaceDst[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            if (ctr_1 > 0)
            {
              fd_edgeFaceDst[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
            }
          }
          {
            fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3];
            if ((1 << level) - 3 > 0)
            {
              fd_edgeFaceDst[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 4];
            }
          }
        }
        {
          fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2];
          if ((1 << level) - 2 > 0)
          {
            fd_edgeFaceDst[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3];
          }
        }
      }
      {
        fd_edgeFaceDst[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1];
        if ((1 << level) - 1 > 0)
        {
          fd_edgeFaceDst[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2];
        }
      }
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + (1 << level) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1];
          if ((1 << level) > 1)
          {
            fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1];
          }
        }
        {
          {
            fd_edgeFaceDst[-(2 / 2) + (1 << level) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2];
            if ((1 << level) > 2)
            {
              fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2];
            }
            fd_edgeFaceDst[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1];
          }
          {
            {
              fd_edgeFaceDst[-(2 / 2) + (1 << level) + 3] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 3] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3];
              if ((1 << level) > 3)
              {
                fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 4] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3];
              }
              fd_edgeFaceDst[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2];
            }
            for (int ctr_1 = 3; ctr_1 < (1 << level) - 4; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (2 / 2) + (1 << level) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1];
              if (ctr_1 + 1 < (1 << level))
              {
                fd_edgeFaceDst[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level)];
              }
            }
            {
              fd_edgeFaceDst[-(2 / 2) + 2*(1 << level) - 3] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 4] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3];
              fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3];
              if ((1 << level) - 4 > 0)
              {
                fd_edgeFaceDst[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 4];
              }
            }
          }
          {
            fd_edgeFaceDst[-(2 / 2) + 2*(1 << level) - 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2];
            fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2];
            if ((1 << level) - 3 > 0)
            {
              fd_edgeFaceDst[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3];
            }
          }
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 2*(1 << level) - 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1];
          fd_edgeFaceDst[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1];
          if ((1 << level) - 2 > 0)
          {
            fd_edgeFaceDst[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2];
          }
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 2*(1 << level) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2];
            if ((1 << level) > 2)
            {
              fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 3] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2];
            }
          }
          {
            {
              fd_edgeFaceDst[-(6 / 2) + 2*(1 << level) + 3] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3];
              if ((1 << level) > 3)
              {
                fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 4] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3];
              }
              fd_edgeFaceDst[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2];
            }
            {
              {
                fd_edgeFaceDst[-(6 / 2) + 2*(1 << level) + 4] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 4] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4];
                if ((1 << level) > 4)
                {
                  fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 5] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 5] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4];
                }
                fd_edgeFaceDst[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 4] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 4] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3];
              }
              for (int ctr_1 = 3; ctr_1 < (1 << level) - 5; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2];
                if (ctr_1 + 2 < (1 << level))
                {
                  fd_edgeFaceDst[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 3] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2];
                }
                if (ctr_1 > 0)
                {
                  fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 1];
                }
              }
              {
                fd_edgeFaceDst[-(6 / 2) + 3*(1 << level) - 3] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 4] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3];
                fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3];
                if ((1 << level) - 5 > 0)
                {
                  fd_edgeFaceDst[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 4];
                }
              }
            }
            {
              fd_edgeFaceDst[-(6 / 2) + 3*(1 << level) - 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2];
              fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2];
              if ((1 << level) - 4 > 0)
              {
                fd_edgeFaceDst[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 3];
              }
            }
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 3*(1 << level) - 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1];
            fd_edgeFaceDst[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1];
            if ((1 << level) - 3 > 0)
            {
              fd_edgeFaceDst[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < (1 << level) - 3; ctr_2 += 1)
        {
          {
            if (ctr_2 > 0)
            {
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            if (ctr_2 < (1 << level))
            {
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
          }
          {
            {
              if (ctr_2 > 0)
              {
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              if (ctr_2 + 1 < (1 << level))
              {
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
            }
            {
              {
                if (ctr_2 > 0)
                {
                  fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
                }
                if (ctr_2 + 2 < (1 << level))
                {
                  fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
                }
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              for (int ctr_1 = 3; ctr_1 < -ctr_2 + (1 << level) - 3; ctr_1 += 1)
              {
                if (ctr_2 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
                }
                if (ctr_1 + ctr_2 < (1 << level))
                {
                  fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
                }
                if (ctr_1 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
                }
              }
              {
                if (ctr_2 > 0)
                {
                  fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 3] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 3] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3];
                }
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 3] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 3] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3];
                if (-ctr_2 + (1 << level) - 3 > 0)
                {
                  fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 4] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 4];
                }
              }
            }
            {
              if (ctr_2 > 0)
              {
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2];
              }
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2];
              if (-ctr_2 + (1 << level) - 2 > 0)
              {
                fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 3] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3];
              }
            }
          }
          {
            if (ctr_2 > 0)
            {
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1];
            }
            fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1];
            if (-ctr_2 + (1 << level) - 1 > 0)
            {
              fd_edgeFaceDst[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 2] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2];
            }
          }
        }
        {
          {
            if ((1 << level) - 3 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
          }
          {
            {
              if ((1 << level) - 3 > 0)
              {
                fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
            }
            {
              {
                if ((1 << level) - 3 > 0)
                {
                  fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
                }
                fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
                fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              for (int ctr_1 = 3; ctr_1 < 0; ctr_1 += 1)
              {
                if ((1 << level) - 3 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
                }
                if (ctr_1 + (1 << level) - 3 < (1 << level))
                {
                  fd_edgeFaceDst[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
                }
                if (ctr_1 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
                }
              }
              {
                if ((1 << level) - 3 > 0)
                {
                  fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
                }
                fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
              }
            }
            {
              if ((1 << level) - 3 > 0)
              {
                fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              }
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
            }
          }
          {
            if ((1 << level) - 3 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 4)*((1 << level) + 1) - (((1 << level) - 4)*((1 << level) - 3) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
            }
            fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
            fd_edgeFaceDst[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
          }
        }
      }
      {
        {
          if ((1 << level) - 2 > 0)
          {
            fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
          }
          fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
        }
        {
          {
            if ((1 << level) - 2 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
            }
            fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
            fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
          }
          {
            {
              if ((1 << level) - 2 > 0)
              {
                fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
              }
              fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
            }
            for (int ctr_1 = 3; ctr_1 < -1; ctr_1 += 1)
            {
              if ((1 << level) - 2 > 0)
              {
                fd_edgeFaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
              }
              if (ctr_1 + (1 << level) - 2 < (1 << level))
              {
                fd_edgeFaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
              }
            }
            {
              if ((1 << level) - 2 > 0)
              {
                fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1];
              }
              fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1];
            }
          }
          {
            if ((1 << level) - 2 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
          }
        }
        {
          if ((1 << level) - 2 > 0)
          {
            fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
          }
          fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
          fd_edgeFaceDst[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
        }
      }
    }
    {
      {
        if ((1 << level) - 1 > 0)
        {
          fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
        }
        fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
      {
        {
          if ((1 << level) - 1 > 0)
          {
            fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1];
          }
          fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)];
        }
        {
          {
            if ((1 << level) - 1 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2];
            }
            fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] = fd_edgeFaceStencil10*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil11*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil12*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeFaceStencil13*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1];
          }
          for (int ctr_1 = 3; ctr_1 < -2; ctr_1 += 1)
          {
            if ((1 << level) - 1 > 0)
            {
              fd_edgeFaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            if (ctr_1 + (1 << level) - 1 < (1 << level))
            {
              fd_edgeFaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
            }
            if (ctr_1 > 0)
            {
              fd_edgeFaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
            }
          }
          {
            if ((1 << level) - 1 > 0)
            {
              fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 2] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 2];
            }
            fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 2] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 2] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 2];
          }
        }
        {
          if ((1 << level) - 1 > 0)
          {
            fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1];
          }
          fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1];
        }
      }
      {
        if ((1 << level) - 1 > 0)
        {
          fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
        }
        fd_edgeFaceDst[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
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


void applyReplace( double*                fd_edgeFaceDst,
                   double*          fd_edgeFaceSrc,
                   double*          fd_edgeFaceStencil,
                   walberla::uint_t level )
{
  apply_2D_macroface_edgedof_to_edgedof_replace(fd_edgeFaceDst, fd_edgeFaceSrc, fd_edgeFaceStencil, static_cast< int64_t >(level));
}

void applyAdd( double* fd_edgeFaceDst, double* fd_edgeFaceSrc, double* fd_edgeFaceStencil, walberla::uint_t level )
{
  for (int ctr_2 = 0; ctr_2 < (1 << level); ctr_2 += 1)
    for (int ctr_1 = 0; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      if (ctr_2 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1]*fd_edgeFaceStencil[2] + fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[1] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[4] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[3] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)]*fd_edgeFaceStencil[0];
      }
      if (ctr_1 + ctr_2 < (1 << level))
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)]*fd_edgeFaceStencil[8] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1]*fd_edgeFaceStencil[7] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[9] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[5] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)]*fd_edgeFaceStencil[6];
      }
      if (ctr_1 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1]*fd_edgeFaceStencil[13] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[10] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1]*fd_edgeFaceStencil[14] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)]*fd_edgeFaceStencil[12] + fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)]*fd_edgeFaceStencil[11];
      }
    }

}


} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
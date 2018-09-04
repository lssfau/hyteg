
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
      fd_edgeFaceDst[-(0 / 2) + 2*(20 / 2) + 3] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 3] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 3] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 3] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 2];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 5] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 5] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (20 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 5];
          fd_edgeFaceDst[-(2 / 2) + (20 / 2) + 5] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 5] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 5];
        }
        for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 5] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 5];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (20 / 2) + 5] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 10] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 5];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(20 / 2) + 5] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 5] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 9] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 4];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 7] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (20 / 2) + 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(20 / 2) + 3] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 7] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7];
          fd_edgeFaceDst[-(2 / 2) + 2*(20 / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6];
        }
      }
      {
        {
          fd_edgeFaceDst[-(6 / 2) + 10] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 10];
          fd_edgeFaceDst[-(6 / 2) + (20 / 2) + 10] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 15] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 10];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (6 / 2) + 10] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 10] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 5] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 6] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(20 / 2) + 10];
        }
        {
          fd_edgeFaceDst[-(6 / 2) + 11] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11];
          fd_edgeFaceDst[-(6 / 2) + 2*(20 / 2) + 11] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (12 / 2) + 15] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 15] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (20 / 2) + 10] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(20 / 2) + 11] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (20 / 2) + 15] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(20 / 2) + 15];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(72 / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 6];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 9] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 9] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (72 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (72 / 2) + 9] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(72 / 2) + 9];
          fd_edgeFaceDst[-(2 / 2) + (72 / 2) + 9] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (72 / 2) + 9] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 9] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(72 / 2) + 10] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 18] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(72 / 2) + 9];
        }
        for (int ctr_1 = 1; ctr_1 < 6; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 9] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 9] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (72 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(72 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (72 / 2) + 9] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(72 / 2) + 9];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (72 / 2) + 9] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (72 / 2) + 9] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 9] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(72 / 2) + 10] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 18] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(72 / 2) + 9];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(72 / 2) + 9] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(72 / 2) + 9] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 9] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (72 / 2) + 9] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 17] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (72 / 2) + 8];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 15] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 15] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (72 / 2) + 6] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(72 / 2) + 7] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (72 / 2) + 15] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(72 / 2) + 15];
          fd_edgeFaceDst[-(2 / 2) + 2*(72 / 2) + 15] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(72 / 2) + 15] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 15] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (72 / 2) + 15] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 23] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (72 / 2) + 14];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 6; ctr_2 += 1)
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
          fd_edgeFaceDst[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] = fd_edgeFaceStencil10*fd_edgeFaceSrc[8*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil11*fd_edgeFaceSrc[8*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil12*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 7] + fd_edgeFaceStencil13*fd_edgeFaceSrc[8*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 15] + fd_edgeFaceStencil14*fd_edgeFaceSrc[8*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 6];
        }
      }
      {
        {
          fd_edgeFaceDst[-(42 / 2) + 54] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(42 / 2) + 54] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(30 / 2) + (72 / 2) + 45] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(30 / 2) + 2*(72 / 2) + 46] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 54];
          fd_edgeFaceDst[-(42 / 2) + (72 / 2) + 54] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(42 / 2) + 54] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 54];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (42 / 2) + 54] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (42 / 2) + 54] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (30 / 2) + (72 / 2) + 45] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (30 / 2) + 2*(72 / 2) + 46] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (42 / 2) + 2*(72 / 2) + 54];
        }
        {
          fd_edgeFaceDst[-(42 / 2) + 55] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(42 / 2) + 55] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(30 / 2) + (72 / 2) + 46] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(30 / 2) + 2*(72 / 2) + 47] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 55] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 55];
          fd_edgeFaceDst[-(42 / 2) + 2*(72 / 2) + 55] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(42 / 2) + 55] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 55] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(56 / 2) + 63] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(42 / 2) + (72 / 2) + 54];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (56 / 2) + 63] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (56 / 2) + 63] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (42 / 2) + (72 / 2) + 54] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (42 / 2) + 2*(72 / 2) + 55] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (56 / 2) + (72 / 2) + 63] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (56 / 2) + 2*(72 / 2) + 63];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(272 / 2) + 15] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 15] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 15] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 15] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 31] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 14];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 17] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 17] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (272 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (272 / 2) + 17] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(272 / 2) + 17];
          fd_edgeFaceDst[-(2 / 2) + (272 / 2) + 17] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (272 / 2) + 17] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 17] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(272 / 2) + 18] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 34] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(272 / 2) + 17];
        }
        for (int ctr_1 = 1; ctr_1 < 14; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 17] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 17] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 17];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (272 / 2) + 17] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 17] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 18] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 34] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 17];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(272 / 2) + 17] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 17] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 17] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 33] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 16];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 31] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 31] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (272 / 2) + 14] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(272 / 2) + 15] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (272 / 2) + 31] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(272 / 2) + 31];
          fd_edgeFaceDst[-(2 / 2) + 2*(272 / 2) + 31] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(272 / 2) + 31] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 31] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (272 / 2) + 31] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 47] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (272 / 2) + 30];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 14; ctr_2 += 1)
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
          fd_edgeFaceDst[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] = fd_edgeFaceStencil10*fd_edgeFaceSrc[16*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil11*fd_edgeFaceSrc[16*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil12*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 15] + fd_edgeFaceStencil13*fd_edgeFaceSrc[16*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 31] + fd_edgeFaceStencil14*fd_edgeFaceSrc[16*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 14];
        }
      }
      {
        {
          fd_edgeFaceDst[-(210 / 2) + 238] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(210 / 2) + 238] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(182 / 2) + (272 / 2) + 221] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(182 / 2) + 2*(272 / 2) + 222] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 238];
          fd_edgeFaceDst[-(210 / 2) + (272 / 2) + 238] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(210 / 2) + 238] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 238];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (210 / 2) + 238] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (210 / 2) + 238] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (182 / 2) + (272 / 2) + 221] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (182 / 2) + 2*(272 / 2) + 222] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (210 / 2) + 2*(272 / 2) + 238];
        }
        {
          fd_edgeFaceDst[-(210 / 2) + 239] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(210 / 2) + 239] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(182 / 2) + (272 / 2) + 222] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(182 / 2) + 2*(272 / 2) + 223] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 239] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 239];
          fd_edgeFaceDst[-(210 / 2) + 2*(272 / 2) + 239] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(210 / 2) + 239] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 239] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(240 / 2) + 255] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(210 / 2) + (272 / 2) + 238];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (240 / 2) + 255] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (240 / 2) + 255] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (210 / 2) + (272 / 2) + 238] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (210 / 2) + 2*(272 / 2) + 239] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (240 / 2) + (272 / 2) + 255] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (240 / 2) + 2*(272 / 2) + 255];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(1056 / 2) + 31] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(1056 / 2) + 31] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 31] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 31] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 63] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (1056 / 2) + 30];
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
        }
        {
          fd_edgeFaceDst[-(930 / 2) + 991] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(930 / 2) + 991] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(870 / 2) + (1056 / 2) + 958] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(870 / 2) + 2*(1056 / 2) + 959] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 991] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991];
          fd_edgeFaceDst[-(930 / 2) + 2*(1056 / 2) + 991] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(930 / 2) + 991] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 991] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(992 / 2) + 1023] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(930 / 2) + (1056 / 2) + 990];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (992 / 2) + 1023] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1023] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 991] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (992 / 2) + (1056 / 2) + 1023] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 2*(1056 / 2) + 1023];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(4160 / 2) + 63] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(4160 / 2) + 63] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 63] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2) + 63] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 127] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (4160 / 2) + 62];
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
        fd_edgeFaceDst[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] = fd_edgeFaceStencil10*fd_edgeFaceSrc[64*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil11*fd_edgeFaceSrc[64*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil12*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 63] + fd_edgeFaceStencil13*fd_edgeFaceSrc[64*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 127] + fd_edgeFaceStencil14*fd_edgeFaceSrc[64*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 62];
      }
    }
    {
      fd_edgeFaceDst[-(4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + 2*(4160 / 2) + 4095];
      }
      fd_edgeFaceDst[-(4032 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4032 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(3906 / 2) + (4160 / 2) + 4030] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(3906 / 2) + 2*(4160 / 2) + 4031] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4032 / 2) + (4160 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4032 / 2) + 2*(4160 / 2) + 4095];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(16512 / 2) + 127] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(16512 / 2) + 127] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 127] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2) + 127] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 255] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (16512 / 2) + 126];
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
        fd_edgeFaceDst[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] = fd_edgeFaceStencil10*fd_edgeFaceSrc[128*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil11*fd_edgeFaceSrc[128*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil12*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 127] + fd_edgeFaceStencil13*fd_edgeFaceSrc[128*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 255] + fd_edgeFaceStencil14*fd_edgeFaceSrc[128*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 126];
      }
    }
    {
      fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (16256 / 2) + 2*(16512 / 2) + 16383];
      }
      fd_edgeFaceDst[-(16256 / 2) + 16383] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16256 / 2) + 16383] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16002 / 2) + (16512 / 2) + 16254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16002 / 2) + 2*(16512 / 2) + 16255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16256 / 2) + (16512 / 2) + 16383] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16256 / 2) + 2*(16512 / 2) + 16383];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(65792 / 2) + 255] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(65792 / 2) + 255] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 255] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 255] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (65792 / 2) + 254];
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
          fd_edgeFaceDst[-(2 / 2) + 2*(65792 / 2) + 511] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(65792 / 2) + 511] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 511] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 511] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (65792 / 2) + 510];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 254; ctr_2 += 1)
      {
        {
          fd_edgeFaceDst[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeFaceStencil2*fd_edgeFaceSrc[257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 256] + fd_edgeFaceStencil3*fd_edgeFaceSrc[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_edgeFaceStencil9*fd_edgeFaceSrc[257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
        }
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 256] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
          fd_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
        }
        {
          fd_edgeFaceDst[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_edgeFaceStencil0*fd_edgeFaceSrc[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_edgeFaceStencil1*fd_edgeFaceSrc[256*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[256*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_edgeFaceStencil4*fd_edgeFaceSrc[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255];
          fd_edgeFaceDst[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] = fd_edgeFaceStencil10*fd_edgeFaceSrc[256*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_edgeFaceStencil11*fd_edgeFaceSrc[256*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_edgeFaceStencil12*fd_edgeFaceSrc[256*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 255] + fd_edgeFaceStencil13*fd_edgeFaceSrc[256*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 511] + fd_edgeFaceStencil14*fd_edgeFaceSrc[256*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 254];
        }
      }
      {
        {
          fd_edgeFaceDst[-(64770 / 2) + 65278] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(64770 / 2) + 65278] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(64262 / 2) + 2*(65792 / 2) + 65022] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65278];
          fd_edgeFaceDst[-(64770 / 2) + (65792 / 2) + 65278] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(64770 / 2) + 65278] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(65280 / 2) + 65535] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65278];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (64770 / 2) + 65278] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 65278] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + (65792 / 2) + 65021] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (64262 / 2) + 2*(65792 / 2) + 65022] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 2*(65792 / 2) + 65278];
        }
        {
          fd_edgeFaceDst[-(64770 / 2) + 65279] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(64770 / 2) + 65279] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(64262 / 2) + (65792 / 2) + 65022] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(64262 / 2) + 2*(65792 / 2) + 65023] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65279] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279];
          fd_edgeFaceDst[-(64770 / 2) + 2*(65792 / 2) + 65279] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(64770 / 2) + 65279] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65279] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(65280 / 2) + 65535] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(64770 / 2) + (65792 / 2) + 65278];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (65280 / 2) + 65535] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 65535] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + (65792 / 2) + 65278] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (64770 / 2) + 2*(65792 / 2) + 65279] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + (65792 / 2) + 65535] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (65280 / 2) + 2*(65792 / 2) + 65535];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(262656 / 2) + 511] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 511] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 511] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 511] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 1023] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 510];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 513] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 513] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 513];
          fd_edgeFaceDst[-(2 / 2) + (262656 / 2) + 513] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 513] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 514] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 1026] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 513];
        }
        for (int ctr_1 = 1; ctr_1 < 510; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 513] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (262656 / 2) + 513] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 514] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1026] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1025] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 512];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 1023] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 1023] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (262656 / 2) + 510] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(262656 / 2) + 511] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 1023] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 1023];
          fd_edgeFaceDst[-(2 / 2) + 2*(262656 / 2) + 1023] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 1023] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 1023] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 1023] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 1535] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 1022];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 1026] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 1026] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 514] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1026];
            fd_edgeFaceDst[-(6 / 2) + (262656 / 2) + 1026] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 1026] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1027] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 1539] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1026];
          }
          for (int ctr_1 = 1; ctr_1 < 509; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 1026] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1026] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 514] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1026];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (262656 / 2) + 1026] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1026] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1027] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1539] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1026];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1026] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1026] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1026] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1538] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1025];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 1535] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 1535] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (262656 / 2) + 1022] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(262656 / 2) + 1023] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1535] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1535];
            fd_edgeFaceDst[-(6 / 2) + 2*(262656 / 2) + 1535] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1535] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 1535] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1535] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 2047] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1534];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 1539] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(12 / 2) + 1539] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1027] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(12 / 2) + 2*(262656 / 2) + 1539];
              fd_edgeFaceDst[-(12 / 2) + (262656 / 2) + 1539] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(12 / 2) + 1539] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(12 / 2) + 2*(262656 / 2) + 1540] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(20 / 2) + 2052] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(12 / 2) + 2*(262656 / 2) + 1539];
            }
            for (int ctr_1 = 1; ctr_1 < 508; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 1539] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1539] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1026] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1027] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1539];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (262656 / 2) + 1539] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1539] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1540] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2052] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1539];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1539] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1539] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1539] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2051] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (262656 / 2) + 1538];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 2047] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(12 / 2) + 2047] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(6 / 2) + (262656 / 2) + 1534] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(6 / 2) + 2*(262656 / 2) + 1535] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(12 / 2) + (262656 / 2) + 2047] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(12 / 2) + 2*(262656 / 2) + 2047];
              fd_edgeFaceDst[-(12 / 2) + 2*(262656 / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(12 / 2) + 2*(262656 / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(12 / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(12 / 2) + (262656 / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(20 / 2) + 2559] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(12 / 2) + (262656 / 2) + 2046];
            }
          }
          {
            for (int ctr_1 = 0; ctr_1 < 508; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (20 / 2) + 2052] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2052] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (262656 / 2) + 1539] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(262656 / 2) + 1540] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (20 / 2) + (262656 / 2) + 2052] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2*(262656 / 2) + 2052];
              if (ctr_1 + 4 < 511)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + (262656 / 2) + 2052] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (20 / 2) + (262656 / 2) + 2052] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2052] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2*(262656 / 2) + 2053] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (30 / 2) + 2565] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2*(262656 / 2) + 2052];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 - (20 / 2) + 2*(262656 / 2) + 2052] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2*(262656 / 2) + 2052] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 2052] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (20 / 2) + (262656 / 2) + 2052] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (30 / 2) + 2564] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (20 / 2) + (262656 / 2) + 2051];
              }
            }
            for (int ctr_2 = 5; ctr_2 < 507; ctr_2 += 1)
              for (int ctr_1 = 0; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
              {
                fd_edgeFaceDst[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 512] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                if (ctr_1 + ctr_2 < 511)
                {
                  fd_edgeFaceDst[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
                }
                if (ctr_1 > 0)
                {
                  fd_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
                }
              }
            for (int ctr_1 = 0; ctr_1 < 5; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (257556 / 2) + 260091] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 260091] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (256542 / 2) + (262656 / 2) + 259578] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (256542 / 2) + 2*(262656 / 2) + 259579] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + (262656 / 2) + 260091] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260091];
              if (ctr_1 + 507 < 511)
              {
                fd_edgeFaceDst[ctr_1 - (257556 / 2) + (262656 / 2) + 260091] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + (262656 / 2) + 260091] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 260091] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260092] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 260604] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260091];
              }
              if (ctr_1 > 0)
              {
                fd_edgeFaceDst[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260091] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260091] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 260091] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + (262656 / 2) + 260091] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 260603] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + (262656 / 2) + 260090];
              }
            }
          }
          {
            {
              fd_edgeFaceDst[-(258572 / 2) + 260604] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(258572 / 2) + 260604] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(257556 / 2) + (262656 / 2) + 260091] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(257556 / 2) + 2*(262656 / 2) + 260092] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260604];
              fd_edgeFaceDst[-(258572 / 2) + (262656 / 2) + 260604] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(258572 / 2) + 260604] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260605] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(259590 / 2) + 261117] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260604];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (258572 / 2) + 260604] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 260604] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + (262656 / 2) + 260091] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (257556 / 2) + 2*(262656 / 2) + 260092] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 2*(262656 / 2) + 260604];
              fd_edgeFaceDst[ctr_1 - (258572 / 2) + (262656 / 2) + 260604] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 260604] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 2*(262656 / 2) + 260605] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + 261117] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 2*(262656 / 2) + 260604];
              fd_edgeFaceDst[ctr_1 - (258572 / 2) + 2*(262656 / 2) + 260604] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 2*(262656 / 2) + 260604] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + 260604] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + 261116] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (258572 / 2) + (262656 / 2) + 260603];
            }
            {
              fd_edgeFaceDst[-(258572 / 2) + 260607] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(258572 / 2) + 260607] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(257556 / 2) + (262656 / 2) + 260094] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(257556 / 2) + 2*(262656 / 2) + 260095] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260607] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260607];
              fd_edgeFaceDst[-(258572 / 2) + 2*(262656 / 2) + 260607] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260607] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(258572 / 2) + 260607] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260607] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(259590 / 2) + 261119] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260606];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(259590 / 2) + 261117] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(259590 / 2) + 261117] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260604] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260605] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261117] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261117];
            fd_edgeFaceDst[-(259590 / 2) + (262656 / 2) + 261117] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261117] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(259590 / 2) + 261117] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261118] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(260610 / 2) + 261630] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261117];
          }
          {
            fd_edgeFaceDst[-(259590 / 2) + 261118] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(259590 / 2) + 261118] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260605] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260606] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261118] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261118];
            fd_edgeFaceDst[-(259590 / 2) + (262656 / 2) + 261118] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261118] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(259590 / 2) + 261118] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261119] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(260610 / 2) + 261631] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261118];
            fd_edgeFaceDst[-(259590 / 2) + 2*(262656 / 2) + 261118] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261118] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(259590 / 2) + 261118] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261118] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(260610 / 2) + 261630] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261117];
          }
          {
            fd_edgeFaceDst[-(259590 / 2) + 261119] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(259590 / 2) + 261119] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(258572 / 2) + (262656 / 2) + 260606] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(258572 / 2) + 2*(262656 / 2) + 260607] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261119] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261119];
            fd_edgeFaceDst[-(259590 / 2) + 2*(262656 / 2) + 261119] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261119] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(259590 / 2) + 261119] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261119] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(260610 / 2) + 261631] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261118];
          }
        }
      }
      {
        {
          fd_edgeFaceDst[-(260610 / 2) + 261630] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(260610 / 2) + 261630] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261117] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261118] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261630];
          fd_edgeFaceDst[-(260610 / 2) + (262656 / 2) + 261630] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(260610 / 2) + 261630] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261630];
        }
        {
          fd_edgeFaceDst[-(260610 / 2) + 261631] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(260610 / 2) + 261631] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(259590 / 2) + (262656 / 2) + 261118] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(259590 / 2) + 2*(262656 / 2) + 261119] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261631] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261631];
          fd_edgeFaceDst[-(260610 / 2) + 2*(262656 / 2) + 261631] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(260610 / 2) + 261631] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261631] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(261632 / 2) + 262143] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(260610 / 2) + (262656 / 2) + 261630];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (261632 / 2) + 262143] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + 262143] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + (262656 / 2) + 261630] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 2*(262656 / 2) + 261631] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + (262656 / 2) + 262143] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + 2*(262656 / 2) + 262143];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(1049600 / 2) + 1023] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(1049600 / 2) + 1023] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 1023] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1023] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (1049600 / 2) + 1022];
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
          fd_edgeFaceDst[-(2 / 2) + 2*(1049600 / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 3071] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2046];
        }
      }
      {
        {
          {
            fd_edgeFaceDst[-(6 / 2) + 2050] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 2050] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 1026] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 2050];
            fd_edgeFaceDst[-(6 / 2) + (1049600 / 2) + 2050] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 2050] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 2051] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(12 / 2) + 3075] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 2050];
          }
          for (int ctr_1 = 1; ctr_1 < 1021; ctr_1 += 1)
          {
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2050] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1026] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2051] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3075] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050];
            fd_edgeFaceDst[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2050] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2050] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3074] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2049];
          }
          {
            fd_edgeFaceDst[-(6 / 2) + 3071] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 3071] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (1049600 / 2) + 2046] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(1049600 / 2) + 2047] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 3071] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 3071];
            fd_edgeFaceDst[-(6 / 2) + 2*(1049600 / 2) + 3071] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 3071] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 3071] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 3071] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(12 / 2) + 4095] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 3070];
          }
        }
        {
          {
            {
              fd_edgeFaceDst[-(12 / 2) + 3075] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(12 / 2) + 3075] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 2051] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(12 / 2) + (1049600 / 2) + 3075] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(12 / 2) + 2*(1049600 / 2) + 3075];
              fd_edgeFaceDst[-(12 / 2) + (1049600 / 2) + 3075] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(12 / 2) + (1049600 / 2) + 3075] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(12 / 2) + 3075] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(12 / 2) + 2*(1049600 / 2) + 3076] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(20 / 2) + 4100] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(12 / 2) + 2*(1049600 / 2) + 3075];
            }
            for (int ctr_1 = 1; ctr_1 < 1020; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 3075] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3075] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1049600 / 2) + 2050] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1049600 / 2) + 2051] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (1049600 / 2) + 3075] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3075];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + (1049600 / 2) + 3075] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (1049600 / 2) + 3075] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3075] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3076] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 4100] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3075];
              fd_edgeFaceDst[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3075] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 2*(1049600 / 2) + 3075] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3075] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (1049600 / 2) + 3075] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (20 / 2) + 4099] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (12 / 2) + (1049600 / 2) + 3074];
            }
            {
              fd_edgeFaceDst[-(12 / 2) + 4095] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(12 / 2) + 4095] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(6 / 2) + (1049600 / 2) + 3070] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(6 / 2) + 2*(1049600 / 2) + 3071] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(12 / 2) + (1049600 / 2) + 4095] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(12 / 2) + 2*(1049600 / 2) + 4095];
              fd_edgeFaceDst[-(12 / 2) + 2*(1049600 / 2) + 4095] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(12 / 2) + 2*(1049600 / 2) + 4095] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(12 / 2) + 4095] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(12 / 2) + (1049600 / 2) + 4095] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(20 / 2) + 5119] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(12 / 2) + (1049600 / 2) + 4094];
            }
          }
          for (int ctr_2 = 4; ctr_2 < 1020; ctr_2 += 1)
          {
            {
              fd_edgeFaceDst[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1024] + fd_edgeFaceStencil3*fd_edgeFaceSrc[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_edgeFaceStencil9*fd_edgeFaceSrc[1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
            }
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1024] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)];
              fd_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1];
            }
            {
              fd_edgeFaceDst[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_edgeFaceStencil0*fd_edgeFaceSrc[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_edgeFaceStencil1*fd_edgeFaceSrc[1024*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2] + fd_edgeFaceStencil2*fd_edgeFaceSrc[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[1024*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_edgeFaceStencil4*fd_edgeFaceSrc[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023];
              fd_edgeFaceDst[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] = fd_edgeFaceStencil10*fd_edgeFaceSrc[1024*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_edgeFaceStencil11*fd_edgeFaceSrc[1024*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_edgeFaceStencil12*fd_edgeFaceSrc[1024*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1023] + fd_edgeFaceStencil13*fd_edgeFaceSrc[1024*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2047] + fd_edgeFaceStencil14*fd_edgeFaceSrc[1024*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 1022];
            }
          }
          {
            {
              fd_edgeFaceDst[-(1041420 / 2) + 1045500] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1041420 / 2) + 1045500] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1039380 / 2) + (1049600 / 2) + 1044475] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1039380 / 2) + 2*(1049600 / 2) + 1044476] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045500];
              fd_edgeFaceDst[-(1041420 / 2) + (1049600 / 2) + 1045500] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(1041420 / 2) + 1045500] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045501] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(1043462 / 2) + 1046525] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045500];
            }
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + 1045500] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 1045500] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (1039380 / 2) + (1049600 / 2) + 1044475] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (1039380 / 2) + 2*(1049600 / 2) + 1044476] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045500];
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 1045500] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045501] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 1046525] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045500];
              fd_edgeFaceDst[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045500] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 2*(1049600 / 2) + 1045500] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + 1045500] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (1043462 / 2) + 1046524] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (1041420 / 2) + (1049600 / 2) + 1045499];
            }
            {
              fd_edgeFaceDst[-(1041420 / 2) + 1045503] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1041420 / 2) + 1045503] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1039380 / 2) + (1049600 / 2) + 1044478] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1039380 / 2) + 2*(1049600 / 2) + 1044479] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045503] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045503];
              fd_edgeFaceDst[-(1041420 / 2) + 2*(1049600 / 2) + 1045503] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045503] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(1041420 / 2) + 1045503] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045503] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(1043462 / 2) + 1046527] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045502];
            }
          }
        }
        {
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046525] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1043462 / 2) + 1046525] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045500] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045501] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046525];
            fd_edgeFaceDst[-(1043462 / 2) + (1049600 / 2) + 1046525] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046525] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(1043462 / 2) + 1046525] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046526] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(1045506 / 2) + 1047550] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046525];
          }
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046526] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1043462 / 2) + 1046526] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045501] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045502] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046526] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046526];
            fd_edgeFaceDst[-(1043462 / 2) + (1049600 / 2) + 1046526] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046526] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(1043462 / 2) + 1046526] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046527] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046526];
            fd_edgeFaceDst[-(1043462 / 2) + 2*(1049600 / 2) + 1046526] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046526] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(1043462 / 2) + 1046526] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046526] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(1045506 / 2) + 1047550] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046525];
          }
          {
            fd_edgeFaceDst[-(1043462 / 2) + 1046527] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(1043462 / 2) + 1046527] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(1041420 / 2) + (1049600 / 2) + 1045502] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(1041420 / 2) + 2*(1049600 / 2) + 1045503] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046527] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046527];
            fd_edgeFaceDst[-(1043462 / 2) + 2*(1049600 / 2) + 1046527] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(1043462 / 2) + 2*(1049600 / 2) + 1046527] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(1043462 / 2) + 1046527] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046527] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(1043462 / 2) + (1049600 / 2) + 1046526];
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
          fd_edgeFaceDst[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(1045506 / 2) + 2*(1049600 / 2) + 1047551] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(1045506 / 2) + 1047551] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047551] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(1047552 / 2) + 1048575] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(1045506 / 2) + (1049600 / 2) + 1047550];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (1047552 / 2) + 1048575] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 1048575] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + (1049600 / 2) + 1047550] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (1045506 / 2) + 2*(1049600 / 2) + 1047551] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + (1049600 / 2) + 1048575] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (1047552 / 2) + 2*(1049600 / 2) + 1048575];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(4196352 / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(4196352 / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 4095] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (4196352 / 2) + 2046];
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
        fd_edgeFaceDst[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] = fd_edgeFaceStencil10*fd_edgeFaceSrc[2048*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil11*fd_edgeFaceSrc[2048*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil12*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2047] + fd_edgeFaceStencil13*fd_edgeFaceSrc[2048*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4095] + fd_edgeFaceStencil14*fd_edgeFaceSrc[2048*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 2046];
      }
    }
    {
      fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + 2*(4196352 / 2) + 4194303];
      }
      fd_edgeFaceDst[-(4192256 / 2) + 4194303] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(4192256 / 2) + 4194303] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(4188162 / 2) + 2*(4196352 / 2) + 4192255] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(4192256 / 2) + (4196352 / 2) + 4194303] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(4192256 / 2) + 2*(4196352 / 2) + 4194303];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(16781312 / 2) + 4095] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 4095] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 4095] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4095] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 8191] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4094];
    }
    {
      {
        {
          fd_edgeFaceDst[-(2 / 2) + 4097] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 4097] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (16781312 / 2) + 4097] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(16781312 / 2) + 4097];
          fd_edgeFaceDst[-(2 / 2) + (16781312 / 2) + 4097] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + (16781312 / 2) + 4097] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 4097] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + 2*(16781312 / 2) + 4098] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 8194] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(2 / 2) + 2*(16781312 / 2) + 4097];
        }
        for (int ctr_1 = 1; ctr_1 < 4094; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 4097] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4098] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 8194] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097];
          fd_edgeFaceDst[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 8193] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4096];
        }
        {
          fd_edgeFaceDst[-(2 / 2) + 8191] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 8191] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (16781312 / 2) + 4094] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(16781312 / 2) + 4095] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (16781312 / 2) + 8191] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + 2*(16781312 / 2) + 8191];
          fd_edgeFaceDst[-(2 / 2) + 2*(16781312 / 2) + 8191] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + 2*(16781312 / 2) + 8191] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 8191] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(2 / 2) + (16781312 / 2) + 8191] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(6 / 2) + 12287] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(2 / 2) + (16781312 / 2) + 8190];
        }
      }
      for (int ctr_2 = 2; ctr_2 < 4094; ctr_2 += 1)
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
          fd_edgeFaceDst[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] = fd_edgeFaceStencil10*fd_edgeFaceSrc[4096*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil11*fd_edgeFaceSrc[4096*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil12*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4095] + fd_edgeFaceStencil13*fd_edgeFaceSrc[4096*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8191] + fd_edgeFaceStencil14*fd_edgeFaceSrc[4096*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 4094];
        }
      }
      {
        {
          fd_edgeFaceDst[-(16764930 / 2) + 16773118] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16764930 / 2) + 16773118] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16756742 / 2) + (16781312 / 2) + 16769021] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16756742 / 2) + 2*(16781312 / 2) + 16769022] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773118];
          fd_edgeFaceDst[-(16764930 / 2) + (16781312 / 2) + 16773118] = fd_edgeFaceStencil5*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil6*fd_edgeFaceSrc[-(16764930 / 2) + 16773118] + fd_edgeFaceStencil7*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil8*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil9*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773118];
        }
        for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
        {
          fd_edgeFaceDst[ctr_1 - (16764930 / 2) + 16773118] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + 16773118] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (16756742 / 2) + (16781312 / 2) + 16769021] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (16756742 / 2) + 2*(16781312 / 2) + 16769022] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + 2*(16781312 / 2) + 16773118];
        }
        {
          fd_edgeFaceDst[-(16764930 / 2) + 16773119] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(16764930 / 2) + 16773119] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(16756742 / 2) + (16781312 / 2) + 16769022] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(16756742 / 2) + 2*(16781312 / 2) + 16769023] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773119] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773119];
          fd_edgeFaceDst[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(16764930 / 2) + 16773119] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773119] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(16773120 / 2) + 16777215] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(16764930 / 2) + (16781312 / 2) + 16773118];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {
      fd_edgeFaceDst[ctr_1 - (16773120 / 2) + 16777215] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + 16777215] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + (16781312 / 2) + 16773118] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (16764930 / 2) + 2*(16781312 / 2) + 16773119] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + (16781312 / 2) + 16777215] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (16773120 / 2) + 2*(16781312 / 2) + 16777215];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(67117056 / 2) + 8191] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(67117056 / 2) + 8191] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 8191] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8191] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 16383] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (67117056 / 2) + 8190];
    }
    for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
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
        fd_edgeFaceDst[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] = fd_edgeFaceStencil10*fd_edgeFaceSrc[8192*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil11*fd_edgeFaceSrc[8192*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil12*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8191] + fd_edgeFaceStencil13*fd_edgeFaceSrc[8192*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16383] + fd_edgeFaceStencil14*fd_edgeFaceSrc[8192*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 8190];
      }
    }
    {
      fd_edgeFaceDst[-(67100672 / 2) + 67108863] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67100672 / 2) + 67108863] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67100672 / 2) + (67117056 / 2) + 67108863] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67100672 / 2) + 2*(67117056 / 2) + 67108863];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (67100672 / 2) + 67108863] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 67108863] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + (67117056 / 2) + 67108863] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 2*(67117056 / 2) + 67108863];
      }
      fd_edgeFaceDst[-(67100672 / 2) + 67108863] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(67100672 / 2) + 67108863] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(67084290 / 2) + 2*(67117056 / 2) + 67100671] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(67100672 / 2) + (67117056 / 2) + 67108863] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(67100672 / 2) + 2*(67117056 / 2) + 67108863];
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
      fd_edgeFaceDst[-(0 / 2) + 2*(268451840 / 2) + 16383] = fd_edgeFaceStencil10*fd_edgeFaceSrc[-(0 / 2) + 2*(268451840 / 2) + 16383] + fd_edgeFaceStencil11*fd_edgeFaceSrc[-(0 / 2) + 16383] + fd_edgeFaceStencil12*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16383] + fd_edgeFaceStencil13*fd_edgeFaceSrc[-(2 / 2) + 32767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[-(0 / 2) + (268451840 / 2) + 16382];
    }
    for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
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
        fd_edgeFaceDst[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] = fd_edgeFaceStencil10*fd_edgeFaceSrc[16384*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil11*fd_edgeFaceSrc[16384*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil12*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16383] + fd_edgeFaceStencil13*fd_edgeFaceSrc[16384*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32767] + fd_edgeFaceStencil14*fd_edgeFaceSrc[16384*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) + 16382];
      }
    }
    {
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(268419072 / 2) + 268435455] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(268419072 / 2) + (268451840 / 2) + 268435455] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(268419072 / 2) + 2*(268451840 / 2) + 268435455];
      for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
      {
        fd_edgeFaceDst[ctr_1 - (268419072 / 2) + 268435455] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 268435455] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + (268451840 / 2) + 268435455] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 2*(268451840 / 2) + 268435455];
      }
      fd_edgeFaceDst[-(268419072 / 2) + 268435455] = fd_edgeFaceStencil0*fd_edgeFaceSrc[-(268419072 / 2) + 268435455] + fd_edgeFaceStencil1*fd_edgeFaceSrc[-(268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeFaceStencil2*fd_edgeFaceSrc[-(268386306 / 2) + 2*(268451840 / 2) + 268419071] + fd_edgeFaceStencil3*fd_edgeFaceSrc[-(268419072 / 2) + (268451840 / 2) + 268435455] + fd_edgeFaceStencil4*fd_edgeFaceSrc[-(268419072 / 2) + 2*(268451840 / 2) + 268435455];
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
  for (int ctr_2 = 0; ctr_2 < (1 << level); ctr_2 += 1)
    for (int ctr_1 = 0; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
    {
      if (ctr_2 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil3*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil4*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
      if (ctr_1 + ctr_2 < (1 << level) - 1)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil5*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeFaceStencil8*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_edgeFaceStencil9*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)];
      }
      if (ctr_1 > 0)
      {
        fd_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] = fd_edgeFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeFaceStencil12*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeFaceStencil13*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_edgeFaceStencil14*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1];
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
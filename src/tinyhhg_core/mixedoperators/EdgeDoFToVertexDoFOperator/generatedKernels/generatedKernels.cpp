
#include "generatedKernels.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_2(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 6] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + (20 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*(20 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_3(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 10] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_4(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 16; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 18] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_5(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 33; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 34] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 32] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 32] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 32] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 65];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 31; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 32)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 68] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 65] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 32] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1056 / 2) + 65] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1056 / 2) + 65] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 34] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 66] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1056 / 2) + 66] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1056 / 2) + 66] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 98];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 30; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 33; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 32)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 34] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 32] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 30 < 32)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (930 / 2) + 1020] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 989] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (870 / 2) + (1056 / 2) + 956] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 989] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 989] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 2*(1056 / 2) + 957] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 957] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (870 / 2) + (1056 / 2) + 957] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 2*(1056 / 2) + 958] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 990] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 990] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1022];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_6(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 65; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 66] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 64] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4160 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4160 / 2) + 64] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4160 / 2) + 64] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4160 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4160 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4160 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 65] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4160 / 2) + 65] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4160 / 2) + 65] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 129];
      }
      {
        for (int ctr_1 = 1; ctr_1 < 62; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (6 / 2) + 132] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 129] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4160 / 2) + 64] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (4160 / 2) + 129] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(4160 / 2) + 129] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4160 / 2) + 65] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 65] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4160 / 2) + 65] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4160 / 2) + 66] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 130] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (4160 / 2) + 130] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(4160 / 2) + 130] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 194];
        }
        for (int ctr_2 = 3; ctr_2 < 62; ctr_2 += 1)
        {
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
          {
            fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 66] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 64] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 64];
          }
        }
        for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (3906 / 2) + 4092] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 4029] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (3782 / 2) + (4160 / 2) + 3964] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + (4160 / 2) + 4029] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 2*(4160 / 2) + 4029] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (3782 / 2) + 2*(4160 / 2) + 3965] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (3782 / 2) + 3965] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (3782 / 2) + (4160 / 2) + 3965] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (3782 / 2) + 2*(4160 / 2) + 3966] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 4030] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + (4160 / 2) + 4030] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (3906 / 2) + 2*(4160 / 2) + 4030] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (4032 / 2) + 4094];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_7(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 130] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 128] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_8(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 256; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 258] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 256] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_9(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 512; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 514] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 512] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_10(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 1025; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 1026] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1024] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1024] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1024] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1049600 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1049600 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 1025] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1049600 / 2) + 1025] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1049600 / 2) + 1025] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2049];
      }
      for (int ctr_2 = 2; ctr_2 < 1023; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1026] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1024] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024];
        }
      }
      for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
      {

      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_11(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 2048; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2050] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2048] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2048];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_12(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 4097; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 4098] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4096] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4096] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4096] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 8193];
      }
      for (int ctr_2 = 2; ctr_2 < 4095; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4098] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4096] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096];
        }
      }
      for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
      {

      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_13(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8194] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8192] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_14(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 16384; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16386] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16384] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < (1 << level) + 1; ctr_1 += 1)
    {

    }
    {
      {
        if ((1 << level) > 1)
        {

        }
        {
          if ((1 << level) > 2)
          {
            fd_p1FaceDst[-(2 / 2) + (1 << level) + 3] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-(0 / 2) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2];
          }
          for (int ctr_1 = 2; ctr_1 < (1 << level) - 2; ctr_1 += 1)
          {
            if (ctr_1 + 1 < (1 << level))
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level)] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 1];
              }
            }
          }
          if ((1 << level) - 2 > 0)
          {
            fd_p1FaceDst[-(2 / 2) + 2*(1 << level)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 3] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-(0 / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-(0 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-(0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1];
          }
        }
      }
      {
        {
          if ((1 << level) > 2)
          {

          }
          {
            if ((1 << level) > 3)
            {
              fd_p1FaceDst[-(6 / 2) + 2*(1 << level) + 5] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 3] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 2*(1 << level) + 3] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 3] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-(12 / 2) + 3*(1 << level) + 3];
            }
            for (int ctr_1 = 2; ctr_1 < (1 << level) - 3; ctr_1 += 1)
            {
              if (ctr_1 + 2 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 2];
                }
              }
            }
            if ((1 << level) - 3 > 0)
            {
              fd_p1FaceDst[-(6 / 2) + 3*(1 << level) + 1] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 2] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 3] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 2] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 2] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) - 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 3*(1 << level) - 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 3*(1 << level) - 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-(12 / 2) + 4*(1 << level) - 1];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < (1 << level) - 2; ctr_2 += 1)
        {
          if (ctr_2 > 0)
          {
            if (ctr_2 < (1 << level))
            {

            }
          }
          {
            if (ctr_2 > 0)
            {
              if (ctr_2 + 1 < (1 << level))
              {
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[(ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[(ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)];
              }
            }
            for (int ctr_1 = 2; ctr_1 < -ctr_2 + (1 << level) - 1; ctr_1 += 1)
            {
              if (ctr_2 > 0)
              {
                if (ctr_1 + ctr_2 < (1 << level))
                {
                  if (ctr_1 > 0)
                  {
                    fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1];
                  }
                }
              }
            }
            if (ctr_2 > 0)
            {
              if (-ctr_2 + (1 << level) - 1 > 0)
              {
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 2] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_2*((1 << level) + 1) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) - 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 2];
              }
            }
          }
          if (ctr_2 > 0)
          {

          }
        }
        {
          if ((1 << level) - 2 > 0)
          {

          }
          {
            if ((1 << level) - 2 > 0)
            {
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)];
            }
            for (int ctr_1 = 2; ctr_1 < 1; ctr_1 += 1)
            {
              if ((1 << level) - 2 > 0)
              {
                if (ctr_1 + (1 << level) - 2 < (1 << level))
                {
                  if (ctr_1 > 0)
                  {
                    fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1];
                  }
                }
              }
            }
            if ((1 << level) - 2 > 0)
            {
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)];
            }
          }
          if ((1 << level) - 2 > 0)
          {

          }
        }
      }
      {
        if ((1 << level) - 1 > 0)
        {

        }
        {
          if ((1 << level) - 1 > 0)
          {

          }
          for (int ctr_1 = 2; ctr_1 < 0; ctr_1 += 1)
          {
            if ((1 << level) - 1 > 0)
            {
              if (ctr_1 + (1 << level) - 1 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1];
                }
              }
            }
          }
          if ((1 << level) - 1 > 0)
          {

          }
        }
        if ((1 << level) - 1 > 0)
        {

        }
      }
    }
    {
      if ((1 << level) > 0)
      {

      }
      {
        if ((1 << level) > 0)
        {

        }
        for (int ctr_1 = 2; ctr_1 < -1; ctr_1 += 1)
        {
          if ((1 << level) > 0)
          {
            if (ctr_1 + (1 << level) < (1 << level))
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + (((1 << level) + 1)*((1 << level) + 1)) - (((1 << level) + 1)*((1 << level) + 2) / 2) - 1];
              }
            }
          }
        }
        if ((1 << level) > 0)
        {

        }
      }
      if ((1 << level) > 0)
      {

      }
    }
  }
}




static void apply_2D_macroface_edgedof_to_vertexdof_replace(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
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



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_2(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 5; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 6] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 4] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 4] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (20 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(20 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 5] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (20 / 2) + 5] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(20 / 2) + 5] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 9] + fd_p1FaceDst[ctr_1 - (2 / 2) + 6];
      }
      fd_p1FaceDst[-(6 / 2) + 13] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[-(6 / 2) + 10] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 5] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 10] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 10] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 6] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[-(2 / 2) + 6] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[-(2 / 2) + (20 / 2) + 6] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[-(2 / 2) + 2*(20 / 2) + 7] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[-(6 / 2) + 11] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[-(6 / 2) + (20 / 2) + 11] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[-(6 / 2) + 2*(20 / 2) + 11] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[-(12 / 2) + 15] + fd_p1FaceDst[-(6 / 2) + 13];
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_3(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 10] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + (72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*(72 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8] + fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_4(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 17; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 18] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 16] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 16] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (272 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(272 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 17] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (272 / 2) + 17] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(272 / 2) + 17] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 33] + fd_p1FaceDst[ctr_1 - (2 / 2) + 18];
      }
      for (int ctr_2 = 2; ctr_2 < 15; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 18] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + (272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*(272 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16] + fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        }
      }
      for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
      {

      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_5(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 33; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 34] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 32] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 32] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 32] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (1056 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(1056 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 65] + fd_p1FaceDst[ctr_1 - (2 / 2) + 34];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 31; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 32)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 68] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 65] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 32] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1056 / 2) + 65] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1056 / 2) + 65] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 33] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 33] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1056 / 2) + 33] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(1056 / 2) + 34] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 66] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (1056 / 2) + 66] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1056 / 2) + 66] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 98] + fd_p1FaceDst[ctr_1 - (6 / 2) + 68];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 30; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 33; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 32)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 34] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 32] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + (1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*(1056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 32] + fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 30 < 32)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (930 / 2) + 1020] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 989] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (870 / 2) + (1056 / 2) + 956] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 989] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 989] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 2*(1056 / 2) + 957] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 957] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (870 / 2) + (1056 / 2) + 957] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (870 / 2) + 2*(1056 / 2) + 958] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 990] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (930 / 2) + (1056 / 2) + 990] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (930 / 2) + 2*(1056 / 2) + 990] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (992 / 2) + 1022] + fd_p1FaceDst[ctr_1 - (930 / 2) + 1020];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_6(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 64; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 66] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 64] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + (4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*(4160 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 64] + fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_7(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 129; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 130] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 128] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16512 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16512 / 2) + 128] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16512 / 2) + 128] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16512 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16512 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16512 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 129] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16512 / 2) + 129] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16512 / 2) + 129] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 257] + fd_p1FaceDst[ctr_1 - (2 / 2) + 130];
      }
      for (int ctr_2 = 2; ctr_2 < 127; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 130] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 128] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + (16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*(16512 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 128] + fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        }
      }
      for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
      {

      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_8(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 256; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 258] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 256] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + (65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*(65792 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 256] + fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_9(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 513; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 514] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 512] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 512] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 512] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (262656 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(262656 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1025] + fd_p1FaceDst[ctr_1 - (2 / 2) + 514];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 511; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 512)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 1028] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1025] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 512] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1025] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1025] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 513] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 513] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (262656 / 2) + 513] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(262656 / 2) + 514] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 1026] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (262656 / 2) + 1026] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(262656 / 2) + 1026] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 1538] + fd_p1FaceDst[ctr_1 - (6 / 2) + 1028];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 510; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 513; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 512)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 514] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 512] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + (262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*(262656 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 512] + fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 510 < 512)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (260610 / 2) + 262140] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 261629] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + (262656 / 2) + 261116] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + (262656 / 2) + 261629] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 2*(262656 / 2) + 261629] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + 2*(262656 / 2) + 261117] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + 261117] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + (262656 / 2) + 261117] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (259590 / 2) + 2*(262656 / 2) + 261118] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 261630] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + (262656 / 2) + 261630] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (260610 / 2) + 2*(262656 / 2) + 261630] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (261632 / 2) + 262142] + fd_p1FaceDst[ctr_1 - (260610 / 2) + 262140];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_10(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1026] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 1024] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + (1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*(1049600 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1024] + fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_11(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 2049; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 2050] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2048] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4196352 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4196352 / 2) + 2048] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4196352 / 2) + 2048] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4196352 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (4196352 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(4196352 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2049] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4196352 / 2) + 2049] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4196352 / 2) + 2049] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 4097] + fd_p1FaceDst[ctr_1 - (2 / 2) + 2050];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 2047; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 2048)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 4100] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 4097] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4196352 / 2) + 2048] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (4196352 / 2) + 4097] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(4196352 / 2) + 4097] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4196352 / 2) + 2049] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2049] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (4196352 / 2) + 2049] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(4196352 / 2) + 2050] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 4098] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (4196352 / 2) + 4098] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(4196352 / 2) + 4098] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 6146] + fd_p1FaceDst[ctr_1 - (6 / 2) + 4100];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 2046; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 2049; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 2048)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2050] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 2048] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + (4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*(4196352 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2048] + fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 2046 < 2048)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (4188162 / 2) + 4194300] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 4192253] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (4184070 / 2) + (4196352 / 2) + 4190204] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + (4196352 / 2) + 4192253] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 2*(4196352 / 2) + 4192253] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (4184070 / 2) + 2*(4196352 / 2) + 4190205] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (4184070 / 2) + 4190205] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (4184070 / 2) + (4196352 / 2) + 4190205] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (4184070 / 2) + 2*(4196352 / 2) + 4190206] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 4192254] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + (4196352 / 2) + 4192254] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (4188162 / 2) + 2*(4196352 / 2) + 4192254] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (4192256 / 2) + 4194302] + fd_p1FaceDst[ctr_1 - (4188162 / 2) + 4194300];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_12(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 4097; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 4098] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4096] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4096] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4096] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (16781312 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(16781312 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 4097] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (16781312 / 2) + 4097] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(16781312 / 2) + 4097] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 8193] + fd_p1FaceDst[ctr_1 - (2 / 2) + 4098];
      }
      for (int ctr_2 = 2; ctr_2 < 4095; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4098] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 4096] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + (16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*(16781312 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4096] + fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
        }
      }
      for (int ctr_1 = 1; ctr_1 < 1; ctr_1 += 1)
      {

      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_13(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 8193; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 8194] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8192] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8192] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8192] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (67117056 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(67117056 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16385] + fd_p1FaceDst[ctr_1 - (2 / 2) + 8194];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 8191; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 8192)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 16388] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16385] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8192] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16385] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16385] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8193] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 8193] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (67117056 / 2) + 8193] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(67117056 / 2) + 8194] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 16386] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (67117056 / 2) + 16386] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(67117056 / 2) + 16386] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 24578] + fd_p1FaceDst[ctr_1 - (6 / 2) + 16388];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 8190; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 8193; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 8192)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8194] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 8192] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + (67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*(67117056 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8192] + fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 8190 < 8192)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (67084290 / 2) + 67108860] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 67100669] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092476] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + (67117056 / 2) + 67100669] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 2*(67117056 / 2) + 67100669] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092477] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 67092477] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + (67117056 / 2) + 67092477] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (67067910 / 2) + 2*(67117056 / 2) + 67092478] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 67100670] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + (67117056 / 2) + 67100670] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (67084290 / 2) + 2*(67117056 / 2) + 67100670] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (67100672 / 2) + 67108862] + fd_p1FaceDst[ctr_1 - (67084290 / 2) + 67108860];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_14(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < 16385; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 16386] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16384] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16384] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16384] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (268451840 / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(268451840 / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 32769] + fd_p1FaceDst[ctr_1 - (2 / 2) + 16386];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 16383; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 16384)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 32772] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 32769] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16384] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (268451840 / 2) + 32769] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(268451840 / 2) + 32769] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16385] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 16385] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (268451840 / 2) + 16385] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(268451840 / 2) + 16386] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 32770] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (268451840 / 2) + 32770] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(268451840 / 2) + 32770] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 49154] + fd_p1FaceDst[ctr_1 - (6 / 2) + 32772];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 16382; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 16385; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 16384)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16386] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 - 1) / 2) - 16384] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + (268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*(268451840 / 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16384] + fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 16382 < 16384)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (268386306 / 2) + 268435452] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 268419069] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + (268451840 / 2) + 268402684] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419069] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419069] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + 2*(268451840 / 2) + 268402685] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + 268402685] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + (268451840 / 2) + 268402685] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (268353542 / 2) + 2*(268451840 / 2) + 268402686] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 268419070] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + (268451840 / 2) + 268419070] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (268386306 / 2) + 2*(268451840 / 2) + 268419070] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (268419072 / 2) + 268435454] + fd_p1FaceDst[ctr_1 - (268386306 / 2) + 268435452];
            }
          }
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_edgedof_to_vertexdof_add_level_any(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
{
  const double fd_edgeToVertexFaceStencil0 = fd_edgeToVertexFaceStencil[0];
  const double fd_edgeToVertexFaceStencil1 = fd_edgeToVertexFaceStencil[1];
  const double fd_edgeToVertexFaceStencil10 = fd_edgeToVertexFaceStencil[10];
  const double fd_edgeToVertexFaceStencil11 = fd_edgeToVertexFaceStencil[11];
  const double fd_edgeToVertexFaceStencil2 = fd_edgeToVertexFaceStencil[2];
  const double fd_edgeToVertexFaceStencil3 = fd_edgeToVertexFaceStencil[3];
  const double fd_edgeToVertexFaceStencil4 = fd_edgeToVertexFaceStencil[4];
  const double fd_edgeToVertexFaceStencil5 = fd_edgeToVertexFaceStencil[5];
  const double fd_edgeToVertexFaceStencil6 = fd_edgeToVertexFaceStencil[6];
  const double fd_edgeToVertexFaceStencil7 = fd_edgeToVertexFaceStencil[7];
  const double fd_edgeToVertexFaceStencil8 = fd_edgeToVertexFaceStencil[8];
  const double fd_edgeToVertexFaceStencil9 = fd_edgeToVertexFaceStencil[9];
  {
    for (int ctr_1 = 0; ctr_1 < (1 << level) + 1; ctr_1 += 1)
    {

    }
    {
      {
        if ((1 << level) > 1)
        {

        }
        for (int ctr_1 = 1; ctr_1 < (1 << level) - 1; ctr_1 += 1)
        {
          if (ctr_1 + 1 < (1 << level))
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level)] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (0 / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (0 / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (0 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 1] + fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2];
            }
          }
        }
      }
      {
        for (int ctr_1 = 0; ctr_1 < (1 << level) - 1; ctr_1 += 1)
        {
          if (ctr_1 + 2 < (1 << level))
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level)] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 - (2 / 2) + (((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 1] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 - (2 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + (1 << level) + 2] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 - (6 / 2) + (((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 - (6 / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 2*(1 << level) + 2] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 2] + fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < (1 << level) - 2; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + (1 << level) + 1; ctr_1 += 1)
          {
            if (ctr_2 > 0)
            {
              if (ctr_1 + ctr_2 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - (ctr_2*(ctr_2 - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - (ctr_2*(ctr_2 + 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)];
                }
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if ((1 << level) - 2 > 0)
          {
            if (ctr_1 + (1 << level) - 2 < (1 << level))
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 1) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 1) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)];
              }
            }
          }
        }
      }
      {
        if ((1 << level) - 1 > 0)
        {

        }
        if ((1 << level) - 1 > 0)
        {

        }
      }
    }
    {
      if ((1 << level) > 0)
      {

      }
      {
        if ((1 << level) > 0)
        {

        }
        for (int ctr_1 = 2; ctr_1 < -1; ctr_1 += 1)
        {
          if ((1 << level) > 0)
          {
            if (ctr_1 + (1 << level) < (1 << level))
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] = fd_edgeToVertexFaceStencil0*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil1*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil10*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - 1] + fd_edgeToVertexFaceStencil11*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) + (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_edgeToVertexFaceStencil2*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil3*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil4*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil5*fd_edgeFaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 1) - (((1 << level) - 1)*(1 << level) / 2) + 2*(((1 << level) + 1)*(1 << level) / 2) + 1] + fd_edgeToVertexFaceStencil6*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil7*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level)] + fd_edgeToVertexFaceStencil8*fd_edgeFaceSrc[ctr_1 + ((1 << level) + 1)*(1 << level) + (((1 << level) + 1)*(1 << level) / 2)] + fd_edgeToVertexFaceStencil9*fd_edgeFaceSrc[ctr_1 + (((1 << level) + 1)*((1 << level) + 1)) - (((1 << level) + 1)*((1 << level) + 2) / 2) - 1] + fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)];
              }
            }
          }
        }
        if ((1 << level) > 0)
        {

        }
      }
      if ((1 << level) > 0)
      {

      }
    }
  }
}




static void apply_2D_macroface_edgedof_to_vertexdof_add(double * fd_edgeFaceSrc, double * fd_edgeToVertexFaceStencil, double * fd_p1FaceDst, int64_t level)
{
  switch( level )
  {
    case 2:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_2(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 3:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_3(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 4:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_4(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 5:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_5(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 6:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_6(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 7:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_7(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 8:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_8(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 9:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_9(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 10:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_10(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 11:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_11(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 12:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_12(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 13:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_13(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    case 14:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_14(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst);
      break;
    default:
      apply_2D_macroface_edgedof_to_vertexdof_add_level_any(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst, level);
      break;
  }
}



void applyFaceReplace( double*                fd_p1FaceDst,
                       double*          fd_edgeFaceSrc,
                       double*          fd_edgeToVertexFaceStencil,
                       walberla::uint_t level )
{
  apply_2D_macroface_edgedof_to_vertexdof_replace(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst, level);
}

void applyFaceAdd( double* fd_p1FaceDst,
                   double* fd_edgeFaceSrc,
                   double* fd_edgeToVertexFaceStencil,
                   walberla::uint_t level )
{
  apply_2D_macroface_edgedof_to_vertexdof_add(fd_edgeFaceSrc, fd_edgeToVertexFaceStencil, fd_p1FaceDst, level);
}


}
}
}

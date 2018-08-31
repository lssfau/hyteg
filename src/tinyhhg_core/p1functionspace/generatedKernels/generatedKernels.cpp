#include "generatedKernels.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_2(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 6] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 5] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 5] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 6];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_3(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 10] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_4(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 16; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 18] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_5(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 32; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 34] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_6(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 64; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 66] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_7(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 130] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_8(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 257; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 258] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 257] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 258] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 259] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 515] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 516];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 255; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 256)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 516] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 258] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 259] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 515] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 516] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 517] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 773] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 774];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 254; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 257; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 256)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 258] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 254 < 256)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (64770 / 2) + 65532] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (64262 / 2) + 65274] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (64262 / 2) + 65275] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65531] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65532] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65533] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (65280 / 2) + 65789] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (65280 / 2) + 65790];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_9(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 512; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 514] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_10(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1026] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_11(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 2049; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 2050] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 2049] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 2050] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 2051] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 4099] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 4100];
      }
      for (int ctr_2 = 2; ctr_2 < 2047; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2050] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_12(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 4097; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 4098] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4097] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4098] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4099] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8195] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8196];
      }
      {
        for (int ctr_1 = 1; ctr_1 < 4094; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (6 / 2) + 8196] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4098] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4099] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8195] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8196] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8197] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 12293] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 12294];
        }
        for (int ctr_2 = 3; ctr_2 < 4094; ctr_2 += 1)
        {
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
          {
            fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4098] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098];
          }
        }
        for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (16764930 / 2) + 16777212] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (16756742 / 2) + 16773114] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (16756742 / 2) + 16773115] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (16764930 / 2) + 16777211] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (16764930 / 2) + 16777212] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (16764930 / 2) + 16777213] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (16773120 / 2) + 16781309] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (16773120 / 2) + 16781310];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_13(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8194] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_14(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 16384; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16386] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_any(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil, int64_t level)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
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
            fd_p1FaceDst[-(2 / 2) + (1 << level) + 3] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 2] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 3] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 4] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 5];
          }
          for (int ctr_1 = 2; ctr_1 < (1 << level) - 2; ctr_1 += 1)
          {
            if (ctr_1 + 1 < (1 << level))
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 2] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 3] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 3] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 4];
              }
            }
          }
          if ((1 << level) - 2 > 0)
          {
            fd_p1FaceDst[-(2 / 2) + 2*(1 << level)] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 2] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 1] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level)] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 1] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 2];
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
              fd_p1FaceDst[-(6 / 2) + 2*(1 << level) + 5] = fd_p1FaceStencil0*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 3] + fd_p1FaceStencil1*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 4] + fd_p1FaceStencil2*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil3*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 5] + fd_p1FaceStencil4*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 6] + fd_p1FaceStencil5*fd_p1FaceSrc[-(12 / 2) + 3*(1 << level) + 6] + fd_p1FaceStencil6*fd_p1FaceSrc[-(12 / 2) + 3*(1 << level) + 7];
            }
            for (int ctr_1 = 2; ctr_1 < (1 << level) - 3; ctr_1 += 1)
            {
              if (ctr_1 + 2 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 2] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 3] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 3] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 5] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 5] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 6];
                }
              }
            }
            if ((1 << level) - 3 > 0)
            {
              fd_p1FaceDst[-(6 / 2) + 3*(1 << level) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level)] + fd_p1FaceStencil2*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level)] + fd_p1FaceStencil3*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[-(12 / 2) + 4*(1 << level) + 2] + fd_p1FaceStencil6*fd_p1FaceSrc[-(12 / 2) + 4*(1 << level) + 3];
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
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[(ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[(ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[(ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[(ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1];
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
                    fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)];
                  }
                }
              }
            }
            if (ctr_2 > 0)
            {
              if (-ctr_2 + (1 << level) - 1 > 0)
              {
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] = fd_p1FaceStencil0*fd_p1FaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + (1 << level) - 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + (1 << level)] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level)] + fd_p1FaceStencil5*fd_p1FaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 2] + fd_p1FaceStencil6*fd_p1FaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 1];
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
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1];
            }
            for (int ctr_1 = 2; ctr_1 < 1; ctr_1 += 1)
            {
              if ((1 << level) - 2 > 0)
              {
                if (ctr_1 + (1 << level) - 2 < (1 << level))
                {
                  if (ctr_1 > 0)
                  {
                    fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)];
                  }
                }
              }
            }
            if ((1 << level) - 2 > 0)
            {
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1];
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
                  fd_p1FaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)];
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
                fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) + 1)*((1 << level) + 2) - (((1 << level) + 1)*((1 << level) + 2) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) + 1)*((1 << level) + 2) - (((1 << level) + 1)*((1 << level) + 2) / 2)];
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




static void apply_2D_macroface_vertexdof_to_vertexdof_replace(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil, int64_t level)
{
  switch( level )
  {
    case 2:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_2(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 3:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_3(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 4:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_4(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 5:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_5(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 6:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_6(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 7:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_7(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 8:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_8(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 9:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_9(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 10:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_10(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 11:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_11(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 12:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_12(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 13:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_13(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 14:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_14(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    default:
      apply_2D_macroface_vertexdof_to_vertexdof_replace_level_any(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil, level);
      break;
  }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_2(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 5; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 6] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 5] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 6] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 7] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 11] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 12] + fd_p1FaceDst[ctr_1 - (2 / 2) + 6];
      }
      for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (6 / 2) + 12] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 6] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 7] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 11] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 12] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 13] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 17] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 18] + fd_p1FaceDst[ctr_1 - (6 / 2) + 12];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_3(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 10] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 9] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 9] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 10] + fd_p1FaceDst[ctr_1 + 10*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_4(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 17; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 18] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 17] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 18] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 19] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 35] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 36] + fd_p1FaceDst[ctr_1 - (2 / 2) + 18];
      }
      {
        for (int ctr_1 = 1; ctr_1 < 14; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (6 / 2) + 36] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 18] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 19] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 35] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 36] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 37] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 53] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 54] + fd_p1FaceDst[ctr_1 - (6 / 2) + 36];
        }
        for (int ctr_2 = 3; ctr_2 < 14; ctr_2 += 1)
        {
          for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
          {
            fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 18] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 17] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 17] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 18] + fd_p1FaceDst[ctr_1 + 18*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
          }
        }
        for (int ctr_1 = 1; ctr_1 < 2; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 - (210 / 2) + 252] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (182 / 2) + 234] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (182 / 2) + 235] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (210 / 2) + 251] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (210 / 2) + 252] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (210 / 2) + 253] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (240 / 2) + 269] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (240 / 2) + 270] + fd_p1FaceDst[ctr_1 - (210 / 2) + 252];
        }
      }
    }
    for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_5(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 32; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 34] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 33] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 33] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 34] + fd_p1FaceDst[ctr_1 + 34*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_6(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 65; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 66] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 65] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 66] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 67] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 131] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 132] + fd_p1FaceDst[ctr_1 - (2 / 2) + 66];
      }
      for (int ctr_2 = 2; ctr_2 < 63; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 66] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 65] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 65] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 66] + fd_p1FaceDst[ctr_1 + 66*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_7(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 130] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 129] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 129] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 130] + fd_p1FaceDst[ctr_1 + 130*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_8(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 257; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 258] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 257] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 258] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 259] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 515] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 516] + fd_p1FaceDst[ctr_1 - (2 / 2) + 258];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 255; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 256)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 516] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 258] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 259] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 515] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 516] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 517] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 773] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 774] + fd_p1FaceDst[ctr_1 - (6 / 2) + 516];
            }
          }
        }
        for (int ctr_2 = 3; ctr_2 < 254; ctr_2 += 1)
          for (int ctr_1 = 0; ctr_1 < -ctr_2 + 257; ctr_1 += 1)
          {
            if (ctr_1 + ctr_2 < 256)
            {
              if (ctr_1 > 0)
              {
                fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 258] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 257] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 257] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 258] + fd_p1FaceDst[ctr_1 + 258*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 254 < 256)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (64770 / 2) + 65532] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (64262 / 2) + 65274] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (64262 / 2) + 65275] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65531] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65532] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (64770 / 2) + 65533] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (65280 / 2) + 65789] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (65280 / 2) + 65790] + fd_p1FaceDst[ctr_1 - (64770 / 2) + 65532];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_9(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 513; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 514] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 513] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 514] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 515] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 1027] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 1028] + fd_p1FaceDst[ctr_1 - (2 / 2) + 514];
      }
      {
        for (int ctr_1 = 0; ctr_1 < 511; ctr_1 += 1)
        {
          if (ctr_1 + 2 < 512)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (6 / 2) + 1028] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + 514] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + 515] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 1027] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 1028] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 1029] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 1541] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 1542] + fd_p1FaceDst[ctr_1 - (6 / 2) + 1028];
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
                fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 514] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 513] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 513] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 514] + fd_p1FaceDst[ctr_1 + 514*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
              }
            }
          }
        for (int ctr_1 = 0; ctr_1 < 3; ctr_1 += 1)
        {
          if (ctr_1 + 510 < 512)
          {
            if (ctr_1 > 0)
            {
              fd_p1FaceDst[ctr_1 - (260610 / 2) + 262140] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (259590 / 2) + 261626] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (259590 / 2) + 261627] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (260610 / 2) + 262139] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (260610 / 2) + 262140] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (260610 / 2) + 262141] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (261632 / 2) + 262653] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (261632 / 2) + 262654] + fd_p1FaceDst[ctr_1 - (260610 / 2) + 262140];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_10(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1026] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 1025] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1025] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1026] + fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_11(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 2048; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2050] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 2049] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2049] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 2050] + fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_12(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 4097; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 4098] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4097] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4098] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 4099] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8195] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 8196] + fd_p1FaceDst[ctr_1 - (2 / 2) + 4098];
      }
      for (int ctr_2 = 2; ctr_2 < 4095; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4098] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 4097] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4097] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 4098] + fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_13(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
    {

    }
    for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
    {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8194] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 8193] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8193] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 8194] + fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
      }
    }
    for (int ctr_1 = 1; ctr_1 < 0; ctr_1 += 1)
    {

    }
  }
}



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_14(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 0; ctr_1 < 16385; ctr_1 += 1)
    {

    }
    {
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
        fd_p1FaceDst[ctr_1 - (2 / 2) + 16386] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + 16385] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + 16386] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + 16387] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 32771] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 32772] + fd_p1FaceDst[ctr_1 - (2 / 2) + 16386];
      }
      for (int ctr_2 = 2; ctr_2 < 16383; ctr_2 += 1)
      {
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
        {
          fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16386] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 - 1) / 2) - 16385] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16385] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 16386] + fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (ctr_2*(ctr_2 + 1) / 2)];
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



static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil, int64_t level)
{
  const double fd_p1FaceStencil0 = fd_p1FaceStencil[0];
  const double fd_p1FaceStencil1 = fd_p1FaceStencil[1];
  const double fd_p1FaceStencil2 = fd_p1FaceStencil[2];
  const double fd_p1FaceStencil3 = fd_p1FaceStencil[3];
  const double fd_p1FaceStencil4 = fd_p1FaceStencil[4];
  const double fd_p1FaceStencil5 = fd_p1FaceStencil[5];
  const double fd_p1FaceStencil6 = fd_p1FaceStencil[6];
  {
    for (int ctr_1 = 1; ctr_1 < (1 << level); ctr_1 += 1)
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
            fd_p1FaceDst[-(2 / 2) + (1 << level) + 3] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 2] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 3] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 4] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 5] + fd_p1FaceDst[-(2 / 2) + (1 << level) + 3];
          }
          {
            if ((1 << level) > 3)
            {
              fd_p1FaceDst[-(2 / 2) + (1 << level) + 4] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + 2] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + 3] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 3] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 4] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 5] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 5] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 6] + fd_p1FaceDst[-(2 / 2) + (1 << level) + 4];
            }
            for (int ctr_1 = 3; ctr_1 < (1 << level) - 3; ctr_1 += 1)
            {
              if (ctr_1 + 1 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (0 / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (0 / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 2] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 3] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 3] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 4] + fd_p1FaceDst[ctr_1 - (2 / 2) + (1 << level) + 2];
                }
              }
            }
            if ((1 << level) - 3 > 0)
            {
              fd_p1FaceDst[-(2 / 2) + 2*(1 << level) - 1] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 3] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 2] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 2] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level)] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level)] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 1] + fd_p1FaceDst[-(2 / 2) + 2*(1 << level) - 1];
            }
          }
          if ((1 << level) - 2 > 0)
          {
            fd_p1FaceDst[-(2 / 2) + 2*(1 << level)] = fd_p1FaceStencil0*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 2] + fd_p1FaceStencil1*fd_p1FaceSrc[-(0 / 2) + (1 << level) - 1] + fd_p1FaceStencil2*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level)] + fd_p1FaceStencil4*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 1] + fd_p1FaceStencil6*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 2] + fd_p1FaceDst[-(2 / 2) + 2*(1 << level)];
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
              fd_p1FaceDst[-(6 / 2) + 2*(1 << level) + 5] = fd_p1FaceStencil0*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 3] + fd_p1FaceStencil1*fd_p1FaceSrc[-(2 / 2) + (1 << level) + 4] + fd_p1FaceStencil2*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil3*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 5] + fd_p1FaceStencil4*fd_p1FaceSrc[-(6 / 2) + 2*(1 << level) + 6] + fd_p1FaceStencil5*fd_p1FaceSrc[-(12 / 2) + 3*(1 << level) + 6] + fd_p1FaceStencil6*fd_p1FaceSrc[-(12 / 2) + 3*(1 << level) + 7] + fd_p1FaceDst[-(6 / 2) + 2*(1 << level) + 5];
            }
            for (int ctr_1 = 2; ctr_1 < (1 << level) - 3; ctr_1 += 1)
            {
              if (ctr_1 + 2 < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 2] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 - (2 / 2) + (1 << level) + 3] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 3] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 4] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 - (6 / 2) + 2*(1 << level) + 5] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 5] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 - (12 / 2) + 3*(1 << level) + 6] + fd_p1FaceDst[ctr_1 - (6 / 2) + 2*(1 << level) + 4];
                }
              }
            }
            if ((1 << level) - 3 > 0)
            {
              fd_p1FaceDst[-(6 / 2) + 3*(1 << level) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level) - 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-(2 / 2) + 2*(1 << level)] + fd_p1FaceStencil2*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level)] + fd_p1FaceStencil3*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[-(6 / 2) + 3*(1 << level) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[-(12 / 2) + 4*(1 << level) + 2] + fd_p1FaceStencil6*fd_p1FaceSrc[-(12 / 2) + 4*(1 << level) + 3] + fd_p1FaceDst[-(6 / 2) + 3*(1 << level) + 1];
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
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[(ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[(ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[(ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[(ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + 1] + fd_p1FaceDst[ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1];
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
                    fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2)] + fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - (ctr_2*(ctr_2 + 1) / 2)];
                  }
                }
              }
            }
            if (ctr_2 > 0)
            {
              if (-ctr_2 + (1 << level) - 1 > 0)
              {
                fd_p1FaceDst[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] = fd_p1FaceStencil0*fd_p1FaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + (1 << level) - 1] + fd_p1FaceStencil1*fd_p1FaceSrc[-ctr_2 + (ctr_2 - 1)*((1 << level) + 2) - (ctr_2*(ctr_2 - 1) / 2) + (1 << level)] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 2] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level)] + fd_p1FaceStencil5*fd_p1FaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 2] + fd_p1FaceStencil6*fd_p1FaceSrc[-ctr_2 + (ctr_2 + 1)*((1 << level) + 2) - ((ctr_2 + 1)*(ctr_2 + 2) / 2) + (1 << level) - 1] + fd_p1FaceDst[ctr_2*((1 << level) + 2) - ctr_2 - (ctr_2*(ctr_2 + 1) / 2) + (1 << level) - 1];
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
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1];
            }
            for (int ctr_1 = 2; ctr_1 < 1; ctr_1 += 1)
            {
              if ((1 << level) - 2 > 0)
              {
                if (ctr_1 + (1 << level) - 2 < (1 << level))
                {
                  if (ctr_1 > 0)
                  {
                    fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceDst[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)];
                  }
                }
              }
            }
            if ((1 << level) - 2 > 0)
            {
              fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] = fd_p1FaceStencil0*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 1] + fd_p1FaceStencil1*fd_p1FaceSrc[((1 << level) - 3)*((1 << level) + 2) - (((1 << level) - 3)*((1 << level) - 2) / 2) + 2] + fd_p1FaceStencil2*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil3*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil4*fd_p1FaceSrc[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 2] + fd_p1FaceStencil5*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil6*fd_p1FaceSrc[((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceDst[((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1];
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
          {
            if ((1 << level) - 1 > 0)
            {

            }
            for (int ctr_1 = 3; ctr_1 < -1; ctr_1 += 1)
            {
              if ((1 << level) - 1 > 0)
              {
                if (ctr_1 + (1 << level) - 1 < (1 << level))
                {
                  if (ctr_1 > 0)
                  {
                    fd_p1FaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 2)*((1 << level) + 2) - (((1 << level) - 2)*((1 << level) - 1) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_p1FaceDst[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)];
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
        {
          if ((1 << level) > 0)
          {

          }
          for (int ctr_1 = 3; ctr_1 < -2; ctr_1 += 1)
          {
            if ((1 << level) > 0)
            {
              if (ctr_1 + (1 << level) < (1 << level))
              {
                if (ctr_1 > 0)
                {
                  fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] = fd_p1FaceStencil0*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2)] + fd_p1FaceStencil1*fd_p1FaceSrc[ctr_1 + ((1 << level) - 1)*((1 << level) + 2) - (((1 << level) - 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil2*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) - 1] + fd_p1FaceStencil3*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)] + fd_p1FaceStencil4*fd_p1FaceSrc[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2) + 1] + fd_p1FaceStencil5*fd_p1FaceSrc[ctr_1 + ((1 << level) + 1)*((1 << level) + 2) - (((1 << level) + 1)*((1 << level) + 2) / 2) - 1] + fd_p1FaceStencil6*fd_p1FaceSrc[ctr_1 + ((1 << level) + 1)*((1 << level) + 2) - (((1 << level) + 1)*((1 << level) + 2) / 2)] + fd_p1FaceDst[ctr_1 + ((1 << level) + 2)*(1 << level) - (((1 << level) + 1)*(1 << level) / 2)];
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
      if ((1 << level) > 0)
      {

      }
    }
  }
}


static void apply_2D_macroface_vertexdof_to_vertexdof_add(double * fd_p1FaceDst, double * fd_p1FaceSrc, double * fd_p1FaceStencil, int64_t level)
{
  switch( level )
  {
    case 2:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_2(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 3:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_3(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 4:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_4(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 5:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_5(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 6:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_6(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 7:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_7(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 8:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_8(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 9:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_9(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 10:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_10(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 11:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_11(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 12:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_12(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 13:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_13(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    case 14:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_14(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil);
      break;
    default:
      apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil, level);
      break;
  }
}



void applyReplace( double*                fd_p1FaceDst,
                   double*          fd_p1FaceSrc,
                   double*          fd_p1FaceStencil,
                   walberla::uint_t level )
{
  apply_2D_macroface_vertexdof_to_vertexdof_replace( fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil, level );
}

void applyAdd( double* fd_p1FaceDst, double* fd_p1FaceSrc, double* fd_p1FaceStencil, walberla::uint_t level )
{
  apply_2D_macroface_vertexdof_to_vertexdof_add( fd_p1FaceDst, fd_p1FaceSrc, fd_p1FaceStencil, level );
}
} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
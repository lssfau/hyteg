
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void assign_2D_macroface_edgedof_1_rhs_function_level_2(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_3(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_4(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_5(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_6(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_7(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_8(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_9(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_10(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_11(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_12(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_13(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_14(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_1_rhs_function_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c, int64_t level)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
      {
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = c*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}


void assign_2D_macroface_edgedof_1_rhs_function(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c, int64_t level)
{
    switch( level )
    {
    case 2:
        assign_2D_macroface_edgedof_1_rhs_function_level_2(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 3:
        assign_2D_macroface_edgedof_1_rhs_function_level_3(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 4:
        assign_2D_macroface_edgedof_1_rhs_function_level_4(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 5:
        assign_2D_macroface_edgedof_1_rhs_function_level_5(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 6:
        assign_2D_macroface_edgedof_1_rhs_function_level_6(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 7:
        assign_2D_macroface_edgedof_1_rhs_function_level_7(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 8:
        assign_2D_macroface_edgedof_1_rhs_function_level_8(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 9:
        assign_2D_macroface_edgedof_1_rhs_function_level_9(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 10:
        assign_2D_macroface_edgedof_1_rhs_function_level_10(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 11:
        assign_2D_macroface_edgedof_1_rhs_function_level_11(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 12:
        assign_2D_macroface_edgedof_1_rhs_function_level_12(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 13:
        assign_2D_macroface_edgedof_1_rhs_function_level_13(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    case 14:
        assign_2D_macroface_edgedof_1_rhs_function_level_14(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c);
        break;
    default:
        assign_2D_macroface_edgedof_1_rhs_function_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, c, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
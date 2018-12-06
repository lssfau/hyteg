
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernels.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_2(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 6] - xi_2*fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_4*fd_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_6*fd_p1FaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + fd_p1FaceRhs[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_3(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 10] - xi_2*fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_4*fd_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_6*fd_p1FaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + fd_p1FaceRhs[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_4(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 18] - xi_2*fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_4*fd_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_6*fd_p1FaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + fd_p1FaceRhs[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_5(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 34] - xi_2*fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_4*fd_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_6*fd_p1FaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + fd_p1FaceRhs[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_6(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 66] - xi_2*fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_4*fd_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_6*fd_p1FaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + fd_p1FaceRhs[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_7(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 130] - xi_2*fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_4*fd_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_6*fd_p1FaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + fd_p1FaceRhs[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_8(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 258] - xi_2*fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_4*fd_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_6*fd_p1FaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + fd_p1FaceRhs[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_9(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 514] - xi_2*fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_4*fd_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_6*fd_p1FaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + fd_p1FaceRhs[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_10(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026] - xi_2*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_4*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_6*fd_p1FaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + fd_p1FaceRhs[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_11(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050] - xi_2*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_4*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_6*fd_p1FaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + fd_p1FaceRhs[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_12(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098] - xi_2*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_4*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_6*fd_p1FaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + fd_p1FaceRhs[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_13(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194] - xi_2*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_4*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_6*fd_p1FaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + fd_p1FaceRhs[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_14(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386] - xi_2*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_4*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_6*fd_p1FaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + fd_p1FaceRhs[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]);
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

static void gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_any(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil, int64_t level)
{
   const double xi_0 = fd_p1FaceStencil[0];
   const double xi_1 = fd_p1FaceStencil[3];
   const double xi_8 = 1 / (xi_1);
   const double xi_2 = fd_p1FaceStencil[2];
   const double xi_3 = fd_p1FaceStencil[1];
   const double xi_4 = fd_p1FaceStencil[4];
   const double xi_5 = fd_p1FaceStencil[6];
   const double xi_6 = fd_p1FaceStencil[5];
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
         fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_8*(-xi_0*fd_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*fd_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << level) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_4*fd_p1FaceDst[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*fd_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_6*fd_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << level) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + fd_p1FaceRhs[ctr_1 + ctr_2*((1 << level) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]);
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


void gaussseidel_2D_macroface_vertexdof_to_vertexdof(double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil, int64_t level)
{
    switch( level )
    {
    case 2:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_2(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 3:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_3(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 4:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_4(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 5:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_5(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 6:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_6(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 7:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_7(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 8:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_8(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 9:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_9(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 10:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_10(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 11:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_11(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 12:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_12(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 13:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_13(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    case 14:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_14(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil);
        break;
    default:
        gaussseidel_2D_macroface_vertexdof_to_vertexdof_level_any(fd_p1FaceDst, fd_p1FaceRhs, fd_p1FaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
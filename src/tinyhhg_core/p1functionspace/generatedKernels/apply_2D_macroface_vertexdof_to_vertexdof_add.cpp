
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsVertexToVertexMacroFace2D.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_2(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 6];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_3(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 10];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_4(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 18];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_5(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 34];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_6(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 66];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_7(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 130];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_8(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 258];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_9(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 514];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_10(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_11(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_12(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_13(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_14(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}

static void apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil, int64_t level)
{
   const double xi_0 = _data_p1FaceStencil[2];
   const double xi_1 = _data_p1FaceStencil[5];
   const double xi_2 = _data_p1FaceStencil[0];
   const double xi_3 = _data_p1FaceStencil[3];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_16 = _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = xi_0*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_10 = xi_1*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_11 = xi_2*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_12 = xi_3*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = xi_4*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_14 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_15 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_9;
      }
   }
}


void apply_2D_macroface_vertexdof_to_vertexdof_add(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_2(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 3:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_3(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 4:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_4(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 5:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_5(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 6:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_6(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 7:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_7(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 8:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_8(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 9:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_9(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 10:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_10(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 11:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_11(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 12:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_12(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 13:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_13(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    case 14:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_14(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil);
        break;
    default:
        apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
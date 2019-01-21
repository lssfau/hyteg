
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_edgedof_to_edgedof_add_level_2(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_12*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_12*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] + xi_9*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_12*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] + xi_12*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + xi_1*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] + _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_3(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_12*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_12*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] + xi_9*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_12*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] + xi_12*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + xi_1*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] + _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_4(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_12*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_12*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] + xi_9*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_12*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] + xi_12*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + xi_1*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] + _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_5(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_12*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_12*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] + xi_9*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_12*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] + xi_12*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + xi_1*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] + _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_6(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_12*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_12*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] + xi_9*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_12*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] + xi_12*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + xi_1*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] + _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_7(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_12*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_12*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] + xi_9*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_12*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] + xi_12*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + xi_1*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] + _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_8(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_12*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_12*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] + xi_9*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_12*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] + xi_12*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + xi_1*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] + _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_9(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_12*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_12*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] + xi_9*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_12*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] + xi_12*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + xi_1*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] + _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_10(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_12*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_12*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] + xi_9*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_12*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] + xi_12*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + xi_1*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] + _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_11(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_12*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_12*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] + xi_9*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_12*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] + xi_12*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + xi_1*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] + _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_12(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_12*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_12*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] + xi_9*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_12*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] + xi_12*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + xi_1*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] + _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_13(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_12*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_12*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] + xi_9*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_12*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] + xi_12*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + xi_1*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] + _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_14(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_12*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_12*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] + xi_9*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_12*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] + xi_12*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + xi_1*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] + _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_any(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil, int64_t level)
{
   const double xi_0 = _data_edgeFaceStencil[1];
   const double xi_1 = _data_edgeFaceStencil[0];
   const double xi_2 = _data_edgeFaceStencil[3];
   const double xi_3 = _data_edgeFaceStencil[4];
   const double xi_4 = _data_edgeFaceStencil[2];
   const double xi_5 = _data_edgeFaceStencil[6];
   const double xi_6 = _data_edgeFaceStencil[5];
   const double xi_7 = _data_edgeFaceStencil[9];
   const double xi_8 = _data_edgeFaceStencil[8];
   const double xi_9 = _data_edgeFaceStencil[7];
   const double xi_10 = _data_edgeFaceStencil[14];
   const double xi_11 = _data_edgeFaceStencil[13];
   const double xi_12 = _data_edgeFaceStencil[11];
   const double xi_13 = _data_edgeFaceStencil[12];
   const double xi_14 = _data_edgeFaceStencil[10];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < (1 << level) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_12*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
      }
      for (int ctr_1 = (1 << level) - 1; ctr_1 < (1 << level); ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_12*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << level) - 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << level) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] = xi_5*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_6*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_7*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_8*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + xi_9*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_12*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + (1 << level) - 1; ctr_1 < -ctr_2 + (1 << level); ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] = xi_10*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2)) - 1] + xi_11*_data_edgeFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << level) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + xi_12*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_13*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_14*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))];
      }
   }
   for (int ctr_2 = (1 << level) - 1; ctr_2 < (1 << level); ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_0*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_1*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + xi_2*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << level) + 1)*(1 << level)) / (2))] + xi_3*_data_edgeFaceSrc[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2))] + xi_4*_data_edgeFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << level) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << level) + 1)*(1 << level)) / (2)) + 1] + _data_edgeFaceDst[ctr_1 + ctr_2*((1 << level) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}


void apply_2D_macroface_edgedof_to_edgedof_add(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeFaceStencil, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_edgedof_to_edgedof_add_level_2(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 3:
        apply_2D_macroface_edgedof_to_edgedof_add_level_3(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 4:
        apply_2D_macroface_edgedof_to_edgedof_add_level_4(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 5:
        apply_2D_macroface_edgedof_to_edgedof_add_level_5(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 6:
        apply_2D_macroface_edgedof_to_edgedof_add_level_6(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 7:
        apply_2D_macroface_edgedof_to_edgedof_add_level_7(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 8:
        apply_2D_macroface_edgedof_to_edgedof_add_level_8(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 9:
        apply_2D_macroface_edgedof_to_edgedof_add_level_9(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 10:
        apply_2D_macroface_edgedof_to_edgedof_add_level_10(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 11:
        apply_2D_macroface_edgedof_to_edgedof_add_level_11(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 12:
        apply_2D_macroface_edgedof_to_edgedof_add_level_12(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 13:
        apply_2D_macroface_edgedof_to_edgedof_add_level_13(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    case 14:
        apply_2D_macroface_edgedof_to_edgedof_add_level_14(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil);
        break;
    default:
        apply_2D_macroface_edgedof_to_edgedof_add_level_any(_data_edgeFaceDst, _data_edgeFaceSrc, _data_edgeFaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
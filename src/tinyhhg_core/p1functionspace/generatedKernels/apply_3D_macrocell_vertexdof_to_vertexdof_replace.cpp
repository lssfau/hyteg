
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsVertexToVertexMacroCell3D.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_2(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 4; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 4; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 5) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 7) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 6)*(-ctr_3 + 7)*(-ctr_3 + 8)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 7) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 6)*(-ctr_3 + 7)*(-ctr_3 + 8)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 5) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 7) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 6)*(-ctr_3 + 7)*(-ctr_3 + 8)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 7) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 6)*(-ctr_3 + 7)*(-ctr_3 + 8)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_3(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 8; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 8; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 9) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 11) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 10)*(-ctr_3 + 11)*(-ctr_3 + 12)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 11) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 10)*(-ctr_3 + 11)*(-ctr_3 + 12)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 9) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 11) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 10)*(-ctr_3 + 11)*(-ctr_3 + 12)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 11) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 10)*(-ctr_3 + 11)*(-ctr_3 + 12)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_4(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 16; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 16; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 17) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 19) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 18)*(-ctr_3 + 19)*(-ctr_3 + 20)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 19) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 18)*(-ctr_3 + 19)*(-ctr_3 + 20)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 17) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 19) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 18)*(-ctr_3 + 19)*(-ctr_3 + 20)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 19) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 18)*(-ctr_3 + 19)*(-ctr_3 + 20)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_5(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 32; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 32; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 32; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 33) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 35) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 34)*(-ctr_3 + 35)*(-ctr_3 + 36)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 35) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 34)*(-ctr_3 + 35)*(-ctr_3 + 36)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 33) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 35) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 34)*(-ctr_3 + 35)*(-ctr_3 + 36)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 35) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 34)*(-ctr_3 + 35)*(-ctr_3 + 36)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_6(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 64; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 64; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 64; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 65) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 67) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 66)*(-ctr_3 + 67)*(-ctr_3 + 68)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 67) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 66)*(-ctr_3 + 67)*(-ctr_3 + 68)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 65) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 67) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 66)*(-ctr_3 + 67)*(-ctr_3 + 68)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 67) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 66)*(-ctr_3 + 67)*(-ctr_3 + 68)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_7(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 128; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 128; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 128; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 131) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 130)*(-ctr_3 + 131)*(-ctr_3 + 132)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 131) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 130)*(-ctr_3 + 131)*(-ctr_3 + 132)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 131) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 130)*(-ctr_3 + 131)*(-ctr_3 + 132)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 131) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 130)*(-ctr_3 + 131)*(-ctr_3 + 132)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_8(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 256; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 256; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 256; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 257) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 259) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 258)*(-ctr_3 + 259)*(-ctr_3 + 260)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 259) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 258)*(-ctr_3 + 259)*(-ctr_3 + 260)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 257) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 259) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 258)*(-ctr_3 + 259)*(-ctr_3 + 260)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 259) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 258)*(-ctr_3 + 259)*(-ctr_3 + 260)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_9(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 512; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 512; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 512; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 515) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 514)*(-ctr_3 + 515)*(-ctr_3 + 516)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 515) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 514)*(-ctr_3 + 515)*(-ctr_3 + 516)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 515) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 514)*(-ctr_3 + 515)*(-ctr_3 + 516)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 515) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 514)*(-ctr_3 + 515)*(-ctr_3 + 516)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_10(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 1024; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 1024; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 1024; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 1027) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 1026)*(-ctr_3 + 1027)*(-ctr_3 + 1028)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 1026) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 1025) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1027) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1026)*(-ctr_3 + 1027)*(-ctr_3 + 1028)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 1027) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 1026)*(-ctr_3 + 1027)*(-ctr_3 + 1028)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 1026) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 1025) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1027) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1026)*(-ctr_3 + 1027)*(-ctr_3 + 1028)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 1026) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1025)*(-ctr_3 + 1026)*(-ctr_3 + 1027)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_11(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 2048; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 2048; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 2048; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 2051) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 2050)*(-ctr_3 + 2051)*(-ctr_3 + 2052)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 2050) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 2049) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2051) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2050)*(-ctr_3 + 2051)*(-ctr_3 + 2052)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 2051) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 2050)*(-ctr_3 + 2051)*(-ctr_3 + 2052)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 2050) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 2049) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2051) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2050)*(-ctr_3 + 2051)*(-ctr_3 + 2052)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 2050) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2049)*(-ctr_3 + 2050)*(-ctr_3 + 2051)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_12(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 4096; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 4096; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4096; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 4099) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 4098)*(-ctr_3 + 4099)*(-ctr_3 + 4100)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 4098) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4097) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4099) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4098)*(-ctr_3 + 4099)*(-ctr_3 + 4100)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 4099) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 4098)*(-ctr_3 + 4099)*(-ctr_3 + 4100)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 4098) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4097) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4099) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4098)*(-ctr_3 + 4099)*(-ctr_3 + 4100)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 4098) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4097)*(-ctr_3 + 4098)*(-ctr_3 + 4099)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_13(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 8192; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 8192; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8192; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 8195) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 8194)*(-ctr_3 + 8195)*(-ctr_3 + 8196)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 8194) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8193) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8195) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8194)*(-ctr_3 + 8195)*(-ctr_3 + 8196)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 8195) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 8194)*(-ctr_3 + 8195)*(-ctr_3 + 8196)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 8194) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8193) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8195) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8194)*(-ctr_3 + 8195)*(-ctr_3 + 8196)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 8194) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8193)*(-ctr_3 + 8194)*(-ctr_3 + 8195)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_14(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < 16384; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 16384; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16384; ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 16387) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 16386)*(-ctr_3 + 16387)*(-ctr_3 + 16388)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 16386) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16385) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16387) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16386)*(-ctr_3 + 16387)*(-ctr_3 + 16388)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 16387) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 16386)*(-ctr_3 + 16387)*(-ctr_3 + 16388)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 16386) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16385) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16387) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16386)*(-ctr_3 + 16387)*(-ctr_3 + 16388)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + 16386) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16385)*(-ctr_3 + 16386)*(-ctr_3 + 16387)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_any(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil, int64_t level)
{
   const double xi_0 = _data_p1CellStencil[2];
   const double xi_1 = _data_p1CellStencil[11];
   const double xi_2 = _data_p1CellStencil[23];
   const double xi_3 = _data_p1CellStencil[5];
   const double xi_4 = _data_p1CellStencil[0];
   const double xi_5 = _data_p1CellStencil[9];
   const double xi_6 = _data_p1CellStencil[21];
   const double xi_7 = _data_p1CellStencil[3];
   const double xi_8 = _data_p1CellStencil[12];
   const double xi_9 = _data_p1CellStencil[24];
   const double xi_10 = _data_p1CellStencil[6];
   const double xi_11 = _data_p1CellStencil[1];
   const double xi_12 = _data_p1CellStencil[10];
   const double xi_13 = _data_p1CellStencil[22];
   const double xi_14 = _data_p1CellStencil[4];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_18 = xi_0*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_19 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_20 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_21 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_22 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) + 1];
            const double xi_23 = xi_13*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) + 1];
            const double xi_24 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}


void apply_3D_macrocell_vertexdof_to_vertexdof_replace(double * _data_p1CellDst, double * _data_p1CellSrc, double * const _data_p1CellStencil, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_2(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 3:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_3(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 4:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_4(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 5:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_5(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 6:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_6(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 7:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_7(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 8:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_8(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 9:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_9(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 10:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_10(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 11:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_11(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 12:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_12(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 13:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_13(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    case 14:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_14(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil);
        break;
    default:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_any(_data_p1CellDst, _data_p1CellSrc, _data_p1CellStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg
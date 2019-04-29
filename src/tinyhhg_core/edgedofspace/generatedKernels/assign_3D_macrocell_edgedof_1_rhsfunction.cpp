
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroCell3D.hpp"

namespace hhg {
namespace edgedof {
namespace macrocell {
namespace generated {

static void assign_3D_macrocell_edgedof_1_rhs_function_level_2(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 3 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 3 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 3 - ctr_2; ctr_1 < 4 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 3 - ctr_3; ctr_2 < 4 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 3; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 3 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 3 - ctr_3; ctr_1 < 4 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 3 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4 - ctr_3) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((3 - ctr_3)*(4 - ctr_3)*(5 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 3; ctr_1 < -ctr_2 - ctr_3 + 4; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 3 - ctr_3; ctr_2 < 4 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 3; ctr_3 < 4; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(5 - ctr_3) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_3(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 7 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 7 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 7 - ctr_2; ctr_1 < 8 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 7 - ctr_3; ctr_2 < 8 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 7; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 7 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 7 - ctr_3; ctr_1 < 8 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 7 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 7; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8 - ctr_3) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((7 - ctr_3)*(8 - ctr_3)*(9 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 7; ctr_1 < -ctr_2 - ctr_3 + 8; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 7 - ctr_3; ctr_2 < 8 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 7; ctr_3 < 8; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(9 - ctr_3) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_4(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 15 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 15 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 15 - ctr_2; ctr_1 < 16 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 15 - ctr_3; ctr_2 < 16 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 15; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 15 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 15 - ctr_3; ctr_1 < 16 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 15 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 15; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16 - ctr_3) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((15 - ctr_3)*(16 - ctr_3)*(17 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 15; ctr_1 < -ctr_2 - ctr_3 + 16; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 15 - ctr_3; ctr_2 < 16 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 15; ctr_3 < 16; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(17 - ctr_3) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_5(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 31 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 31 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 31 - ctr_2; ctr_1 < 32 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 31 - ctr_3; ctr_2 < 32 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 31; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 31 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 31 - ctr_3; ctr_1 < 32 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 31 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 31; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(32 - ctr_3) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((31 - ctr_3)*(32 - ctr_3)*(33 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 31; ctr_1 < -ctr_2 - ctr_3 + 32; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 31 - ctr_3; ctr_2 < 32 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 31; ctr_3 < 32; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(33 - ctr_3) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_6(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 63 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 63 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 63 - ctr_2; ctr_1 < 64 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 63 - ctr_3; ctr_2 < 64 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 63; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 63 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 63 - ctr_3; ctr_1 < 64 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 63 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 63; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(64 - ctr_3) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((63 - ctr_3)*(64 - ctr_3)*(65 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 63; ctr_1 < -ctr_2 - ctr_3 + 64; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 63 - ctr_3; ctr_2 < 64 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 63; ctr_3 < 64; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(65 - ctr_3) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_7(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 127 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 127 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 127 - ctr_2; ctr_1 < 128 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 127 - ctr_3; ctr_2 < 128 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 127; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 127 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 127 - ctr_3; ctr_1 < 128 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 127 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 127; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(128 - ctr_3) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((127 - ctr_3)*(128 - ctr_3)*(129 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 127; ctr_1 < -ctr_2 - ctr_3 + 128; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 127 - ctr_3; ctr_2 < 128 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 127; ctr_3 < 128; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(129 - ctr_3) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_8(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 255 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 255 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 255 - ctr_2; ctr_1 < 256 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 255 - ctr_3; ctr_2 < 256 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 255; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 255 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 255 - ctr_3; ctr_1 < 256 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 255 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 255; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(256 - ctr_3) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((255 - ctr_3)*(256 - ctr_3)*(257 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 255; ctr_1 < -ctr_2 - ctr_3 + 256; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 255 - ctr_3; ctr_2 < 256 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 255; ctr_3 < 256; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(257 - ctr_3) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_9(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 511 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 511 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 511 - ctr_2; ctr_1 < 512 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 511 - ctr_3; ctr_2 < 512 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 511; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 511 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 511 - ctr_3; ctr_1 < 512 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 511 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 511; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(512 - ctr_3) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((511 - ctr_3)*(512 - ctr_3)*(513 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 511; ctr_1 < -ctr_2 - ctr_3 + 512; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 511 - ctr_3; ctr_2 < 512 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 511; ctr_3 < 512; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(513 - ctr_3) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_10(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 1023 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 1023 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 1023 - ctr_2; ctr_1 < 1024 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1023 - ctr_3; ctr_2 < 1024 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 1023; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 1023 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 1023 - ctr_3; ctr_1 < 1024 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 1023 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 1023; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(1024 - ctr_3) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1023 - ctr_3)*(1024 - ctr_3)*(1025 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 1023; ctr_1 < -ctr_2 - ctr_3 + 1024; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1023 - ctr_3; ctr_2 < 1024 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1023; ctr_3 < 1024; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(1025 - ctr_3) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_11(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 2047 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 2047 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 2047 - ctr_2; ctr_1 < 2048 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 2047 - ctr_3; ctr_2 < 2048 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 2047; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 2047 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 2047 - ctr_3; ctr_1 < 2048 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 2047 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 2047; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(2048 - ctr_3) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2047 - ctr_3)*(2048 - ctr_3)*(2049 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 2047; ctr_1 < -ctr_2 - ctr_3 + 2048; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 2047 - ctr_3; ctr_2 < 2048 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 2047; ctr_3 < 2048; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(2049 - ctr_3) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_12(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 4095 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 4095 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 4095 - ctr_2; ctr_1 < 4096 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 4095 - ctr_3; ctr_2 < 4096 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 4095; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 4095 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 4095 - ctr_3; ctr_1 < 4096 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 4095 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4095; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(4096 - ctr_3) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4095 - ctr_3)*(4096 - ctr_3)*(4097 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 4095; ctr_1 < -ctr_2 - ctr_3 + 4096; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 4095 - ctr_3; ctr_2 < 4096 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 4095; ctr_3 < 4096; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(4097 - ctr_3) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_13(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 8191 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 8191 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 8191 - ctr_2; ctr_1 < 8192 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 8191 - ctr_3; ctr_2 < 8192 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 8191; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 8191 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 8191 - ctr_3; ctr_1 < 8192 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 8191 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8191; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(8192 - ctr_3) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8191 - ctr_3)*(8192 - ctr_3)*(8193 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 8191; ctr_1 < -ctr_2 - ctr_3 + 8192; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 8191 - ctr_3; ctr_2 < 8192 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 8191; ctr_3 < 8192; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(8193 - ctr_3) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_14(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // vertex 1
            for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 16383 - ctr_3; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < 16383 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // edge 2
            for (int ctr_1 = 16383 - ctr_2; ctr_1 < 16384 - ctr_2; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 16383 - ctr_3; ctr_2 < 16384 - ctr_3; ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < 16383; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < 16383 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // edge 4
            for (int ctr_1 = 16383 - ctr_3; ctr_1 < 16384 - ctr_3; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < 16383 - ctr_3; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16383; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(16384 - ctr_3) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16383 - ctr_3)*(16384 - ctr_3)*(16385 - ctr_3)) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + 16383; ctr_1 < -ctr_2 - ctr_3 + 16384; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
         for (int ctr_2 = 16383 - ctr_3; ctr_2 < 16384 - ctr_3; ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
      }
      for (int ctr_3 = 16383; ctr_3 < 16384; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(16385 - ctr_3) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            }
         }
      }
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_any(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c, int64_t level)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // vertex 1
            for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)) - 1; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // edge 2
            for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
         for (int ctr_2 = -ctr_3 + (1 << (level)) - 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
         {
            // vertex 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < (1 << (level)) - 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < -ctr_3 + (1 << (level)) - 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // edge 4
            for (int ctr_1 = -ctr_3 + (1 << (level)) - 1; ctr_1 < -ctr_3 + (1 << (level)); ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
         for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)) - 1; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)) - 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + (1 << (level)) - 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
         for (int ctr_2 = -ctr_3 + (1 << (level)) - 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
      }
      for (int ctr_3 = (1 << (level)) - 1; ctr_3 < (1 << (level)); ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
               _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            }
         }
      }
   }
}


void assign_3D_macrocell_edgedof_1_rhs_function(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double * RESTRICT _data_edgeCellSrc_X, double * RESTRICT _data_edgeCellSrc_XY, double * RESTRICT _data_edgeCellSrc_XYZ, double * RESTRICT _data_edgeCellSrc_XZ, double * RESTRICT _data_edgeCellSrc_Y, double * RESTRICT _data_edgeCellSrc_YZ, double * RESTRICT _data_edgeCellSrc_Z, double c, int64_t level)
{
    switch( level )
    {
    case 2:
        assign_3D_macrocell_edgedof_1_rhs_function_level_2(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 3:
        assign_3D_macrocell_edgedof_1_rhs_function_level_3(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 4:
        assign_3D_macrocell_edgedof_1_rhs_function_level_4(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 5:
        assign_3D_macrocell_edgedof_1_rhs_function_level_5(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 6:
        assign_3D_macrocell_edgedof_1_rhs_function_level_6(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 7:
        assign_3D_macrocell_edgedof_1_rhs_function_level_7(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 8:
        assign_3D_macrocell_edgedof_1_rhs_function_level_8(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 9:
        assign_3D_macrocell_edgedof_1_rhs_function_level_9(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 10:
        assign_3D_macrocell_edgedof_1_rhs_function_level_10(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 11:
        assign_3D_macrocell_edgedof_1_rhs_function_level_11(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 12:
        assign_3D_macrocell_edgedof_1_rhs_function_level_12(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 13:
        assign_3D_macrocell_edgedof_1_rhs_function_level_13(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    case 14:
        assign_3D_macrocell_edgedof_1_rhs_function_level_14(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c);
        break;
    default:
        assign_3D_macrocell_edgedof_1_rhs_function_level_any(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, c, level);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace edgedof
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroCell3D.hpp"

namespace hhg {
namespace edgedof {
namespace macrocell {
namespace generated {

static void assign_3D_macrocell_edgedof_1_rhs_function_level_2(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 3; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 3; ctr_2 < -ctr_3 + 4; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 3; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 3; ctr_1 < -ctr_3 + 4; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 3; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 3; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 3; ctr_1 < -ctr_2 - ctr_3 + 4; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 3; ctr_2 < -ctr_3 + 4; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_3(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 7; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 7; ctr_2 < -ctr_3 + 8; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 7; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 7; ctr_1 < -ctr_3 + 8; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 7; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 7; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 7; ctr_1 < -ctr_2 - ctr_3 + 8; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 7; ctr_2 < -ctr_3 + 8; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_4(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 15; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 15; ctr_2 < -ctr_3 + 16; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 15; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 15; ctr_1 < -ctr_3 + 16; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 15; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 15; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 15; ctr_1 < -ctr_2 - ctr_3 + 16; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 15; ctr_2 < -ctr_3 + 16; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_5(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 31; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 31; ctr_2 < -ctr_3 + 32; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 31; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 31; ctr_1 < -ctr_3 + 32; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 31; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 31; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 31; ctr_1 < -ctr_2 - ctr_3 + 32; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 31; ctr_2 < -ctr_3 + 32; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_6(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 63; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 63; ctr_2 < -ctr_3 + 64; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 63; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 63; ctr_1 < -ctr_3 + 64; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 63; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 63; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 63; ctr_1 < -ctr_2 - ctr_3 + 64; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 63; ctr_2 < -ctr_3 + 64; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_7(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 127; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 127; ctr_2 < -ctr_3 + 128; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 127; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 127; ctr_1 < -ctr_3 + 128; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 127; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 127; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 127; ctr_1 < -ctr_2 - ctr_3 + 128; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 127; ctr_2 < -ctr_3 + 128; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_8(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 255; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 255; ctr_2 < -ctr_3 + 256; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 255; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 255; ctr_1 < -ctr_3 + 256; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 255; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 255; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 255; ctr_1 < -ctr_2 - ctr_3 + 256; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 255; ctr_2 < -ctr_3 + 256; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_9(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 511; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 511; ctr_2 < -ctr_3 + 512; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 511; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 511; ctr_1 < -ctr_3 + 512; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 511; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 511; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 511; ctr_1 < -ctr_2 - ctr_3 + 512; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 511; ctr_2 < -ctr_3 + 512; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_10(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 1023; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 1023; ctr_2 < -ctr_3 + 1024; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 1023; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 1023; ctr_1 < -ctr_3 + 1024; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 1023; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 1023; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 1024) + ((1073740800) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1023)*(-ctr_3 + 1024)*(-ctr_3 + 1025)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 1023; ctr_1 < -ctr_2 - ctr_3 + 1024; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 1023; ctr_2 < -ctr_3 + 1024; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 1025) + ((1076889600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 1024)*(-ctr_3 + 1025)*(-ctr_3 + 1026)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_11(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 2047; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 2047; ctr_2 < -ctr_3 + 2048; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 2047; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 2047; ctr_1 < -ctr_3 + 2048; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 2047; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 2047; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 2048) + ((8589932544) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2047)*(-ctr_3 + 2048)*(-ctr_3 + 2049)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 2047; ctr_1 < -ctr_2 - ctr_3 + 2048; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 2047; ctr_2 < -ctr_3 + 2048; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 2049) + ((8602521600) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 2048)*(-ctr_3 + 2049)*(-ctr_3 + 2050)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_12(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 4095; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 4095; ctr_2 < -ctr_3 + 4096; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 4095; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 4095; ctr_1 < -ctr_3 + 4096; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 4095; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4095; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4096) + ((68719472640) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4095)*(-ctr_3 + 4096)*(-ctr_3 + 4097)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 4095; ctr_1 < -ctr_2 - ctr_3 + 4096; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 4095; ctr_2 < -ctr_3 + 4096; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 4097) + ((68769816576) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4096)*(-ctr_3 + 4097)*(-ctr_3 + 4098)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_13(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 8191; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 8191; ctr_2 < -ctr_3 + 8192; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 8191; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 8191; ctr_1 < -ctr_3 + 8192; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 8191; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8191; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8192) + ((549755805696) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8191)*(-ctr_3 + 8192)*(-ctr_3 + 8193)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 8191; ctr_1 < -ctr_2 - ctr_3 + 8192; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 8191; ctr_2 < -ctr_3 + 8192; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 8193) + ((549957156864) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8192)*(-ctr_3 + 8193)*(-ctr_3 + 8194)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_14(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c)
{
   for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // vertex 0
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // edge 0
         for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // vertex 1
         for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 16383; ctr_2 += 1)
      {
         // edge 1
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // face 0
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // edge 2
         for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 16383; ctr_2 < -ctr_3 + 16384; ctr_2 += 1)
      {
         // vertex 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + 16383; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // edge 4
         for (int ctr_1 = -ctr_3 + 16383; ctr_1 < -ctr_3 + 16384; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 16383; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16383; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))] = c*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16384) + ((4398046494720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16383)*(-ctr_3 + 16384)*(-ctr_3 + 16385)) / (6))];
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + 16383; ctr_1 < -ctr_2 - ctr_3 + 16384; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
         }
      }
      for (int ctr_2 = -ctr_3 + 16383; ctr_2 < -ctr_3 + 16384; ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
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
            _data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
            _data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))] = c*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 16385) + ((4398851850240) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16384)*(-ctr_3 + 16385)*(-ctr_3 + 16386)) / (6))];
         }
      }
   }
   {
      
   }
}

static void assign_3D_macrocell_edgedof_1_rhs_function_level_any(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c, int64_t level)
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
   {
      
   }
}


void assign_3D_macrocell_edgedof_1_rhs_function(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c, int64_t level)
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsVertexToVertexMacroCell3D.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_2(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 4; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 4 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(6 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(5 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (7 - ctr_3)*(ctr_2 + 1) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((6 - ctr_3)*(7 - ctr_3)*(8 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (6 - ctr_3)*(ctr_2 + 1) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (6 - ctr_3)*(ctr_2 - 1) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (5 - ctr_3)*(ctr_2 - 1) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(7 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((6 - ctr_3)*(7 - ctr_3)*(8 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(6 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(5 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (7 - ctr_3)*(ctr_2 + 1) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((6 - ctr_3)*(7 - ctr_3)*(8 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (6 - ctr_3)*(ctr_2 + 1) + ((210) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (6 - ctr_3)*(ctr_2 - 1) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (5 - ctr_3)*(ctr_2 - 1) + ((210) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4 - ctr_3)*(5 - ctr_3)*(6 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(7 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((6 - ctr_3)*(7 - ctr_3)*(8 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(6 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(6 - ctr_3) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((5 - ctr_3)*(6 - ctr_3)*(7 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_3(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 8; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 8 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(10 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(9 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (11 - ctr_3)*(ctr_2 + 1) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((10 - ctr_3)*(11 - ctr_3)*(12 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (10 - ctr_3)*(ctr_2 + 1) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (10 - ctr_3)*(ctr_2 - 1) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (9 - ctr_3)*(ctr_2 - 1) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(11 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((10 - ctr_3)*(11 - ctr_3)*(12 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(10 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(9 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (11 - ctr_3)*(ctr_2 + 1) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((10 - ctr_3)*(11 - ctr_3)*(12 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (10 - ctr_3)*(ctr_2 + 1) + ((990) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (10 - ctr_3)*(ctr_2 - 1) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (9 - ctr_3)*(ctr_2 - 1) + ((990) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8 - ctr_3)*(9 - ctr_3)*(10 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(11 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((10 - ctr_3)*(11 - ctr_3)*(12 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(10 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(10 - ctr_3) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((9 - ctr_3)*(10 - ctr_3)*(11 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_4(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 16; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 16 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(18 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(17 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (19 - ctr_3)*(ctr_2 + 1) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((18 - ctr_3)*(19 - ctr_3)*(20 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (18 - ctr_3)*(ctr_2 + 1) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (18 - ctr_3)*(ctr_2 - 1) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (17 - ctr_3)*(ctr_2 - 1) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(19 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((18 - ctr_3)*(19 - ctr_3)*(20 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(18 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(17 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (19 - ctr_3)*(ctr_2 + 1) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((18 - ctr_3)*(19 - ctr_3)*(20 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (18 - ctr_3)*(ctr_2 + 1) + ((5814) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (18 - ctr_3)*(ctr_2 - 1) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (17 - ctr_3)*(ctr_2 - 1) + ((5814) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16 - ctr_3)*(17 - ctr_3)*(18 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(19 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((18 - ctr_3)*(19 - ctr_3)*(20 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(18 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(18 - ctr_3) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((17 - ctr_3)*(18 - ctr_3)*(19 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_5(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 32; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 32 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 32; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(34 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(33 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (35 - ctr_3)*(ctr_2 + 1) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((34 - ctr_3)*(35 - ctr_3)*(36 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (34 - ctr_3)*(ctr_2 + 1) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (34 - ctr_3)*(ctr_2 - 1) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (33 - ctr_3)*(ctr_2 - 1) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(35 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((34 - ctr_3)*(35 - ctr_3)*(36 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(34 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(33 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (35 - ctr_3)*(ctr_2 + 1) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((34 - ctr_3)*(35 - ctr_3)*(36 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (34 - ctr_3)*(ctr_2 + 1) + ((39270) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (34 - ctr_3)*(ctr_2 - 1) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (33 - ctr_3)*(ctr_2 - 1) + ((39270) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((32 - ctr_3)*(33 - ctr_3)*(34 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(35 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((34 - ctr_3)*(35 - ctr_3)*(36 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(34 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(34 - ctr_3) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((33 - ctr_3)*(34 - ctr_3)*(35 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_6(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 64; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 64 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 64; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(66 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(65 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (67 - ctr_3)*(ctr_2 + 1) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((66 - ctr_3)*(67 - ctr_3)*(68 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (66 - ctr_3)*(ctr_2 + 1) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (66 - ctr_3)*(ctr_2 - 1) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (65 - ctr_3)*(ctr_2 - 1) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(67 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((66 - ctr_3)*(67 - ctr_3)*(68 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(66 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(65 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (67 - ctr_3)*(ctr_2 + 1) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((66 - ctr_3)*(67 - ctr_3)*(68 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (66 - ctr_3)*(ctr_2 + 1) + ((287430) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (66 - ctr_3)*(ctr_2 - 1) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (65 - ctr_3)*(ctr_2 - 1) + ((287430) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((64 - ctr_3)*(65 - ctr_3)*(66 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(67 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((66 - ctr_3)*(67 - ctr_3)*(68 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(66 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(66 - ctr_3) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((65 - ctr_3)*(66 - ctr_3)*(67 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_7(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 128; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 128 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 128; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(130 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(129 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (131 - ctr_3)*(ctr_2 + 1) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((130 - ctr_3)*(131 - ctr_3)*(132 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (130 - ctr_3)*(ctr_2 + 1) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (130 - ctr_3)*(ctr_2 - 1) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (129 - ctr_3)*(ctr_2 - 1) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(131 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((130 - ctr_3)*(131 - ctr_3)*(132 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(130 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(129 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (131 - ctr_3)*(ctr_2 + 1) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((130 - ctr_3)*(131 - ctr_3)*(132 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (130 - ctr_3)*(ctr_2 + 1) + ((2196870) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (130 - ctr_3)*(ctr_2 - 1) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (129 - ctr_3)*(ctr_2 - 1) + ((2196870) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((128 - ctr_3)*(129 - ctr_3)*(130 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(131 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((130 - ctr_3)*(131 - ctr_3)*(132 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(130 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(130 - ctr_3) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((129 - ctr_3)*(130 - ctr_3)*(131 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_8(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 256; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 256 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 256; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(258 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(257 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (259 - ctr_3)*(ctr_2 + 1) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((258 - ctr_3)*(259 - ctr_3)*(260 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (258 - ctr_3)*(ctr_2 + 1) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (258 - ctr_3)*(ctr_2 - 1) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (257 - ctr_3)*(ctr_2 - 1) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(259 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((258 - ctr_3)*(259 - ctr_3)*(260 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(258 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(257 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (259 - ctr_3)*(ctr_2 + 1) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((258 - ctr_3)*(259 - ctr_3)*(260 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (258 - ctr_3)*(ctr_2 + 1) + ((17173254) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (258 - ctr_3)*(ctr_2 - 1) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (257 - ctr_3)*(ctr_2 - 1) + ((17173254) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((256 - ctr_3)*(257 - ctr_3)*(258 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(259 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((258 - ctr_3)*(259 - ctr_3)*(260 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(258 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(258 - ctr_3) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((257 - ctr_3)*(258 - ctr_3)*(259 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_9(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 512; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 512 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 512; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(514 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(513 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (515 - ctr_3)*(ctr_2 + 1) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((514 - ctr_3)*(515 - ctr_3)*(516 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (514 - ctr_3)*(ctr_2 + 1) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (514 - ctr_3)*(ctr_2 - 1) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (513 - ctr_3)*(ctr_2 - 1) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(515 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((514 - ctr_3)*(515 - ctr_3)*(516 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(514 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(513 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (515 - ctr_3)*(ctr_2 + 1) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((514 - ctr_3)*(515 - ctr_3)*(516 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (514 - ctr_3)*(ctr_2 + 1) + ((135796230) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (514 - ctr_3)*(ctr_2 - 1) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (513 - ctr_3)*(ctr_2 - 1) + ((135796230) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((512 - ctr_3)*(513 - ctr_3)*(514 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(515 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((514 - ctr_3)*(515 - ctr_3)*(516 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(514 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(514 - ctr_3) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((513 - ctr_3)*(514 - ctr_3)*(515 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_10(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 1024; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 1024 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 1024; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(1026 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(1025 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (1027 - ctr_3)*(ctr_2 + 1) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((1026 - ctr_3)*(1027 - ctr_3)*(1028 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (1026 - ctr_3)*(ctr_2 + 1) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (1026 - ctr_3)*(ctr_2 - 1) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (1025 - ctr_3)*(ctr_2 - 1) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(1027 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1026 - ctr_3)*(1027 - ctr_3)*(1028 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(1026 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(1025 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (1027 - ctr_3)*(ctr_2 + 1) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((1026 - ctr_3)*(1027 - ctr_3)*(1028 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (1026 - ctr_3)*(ctr_2 + 1) + ((1080044550) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (1026 - ctr_3)*(ctr_2 - 1) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (1025 - ctr_3)*(ctr_2 - 1) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((1024 - ctr_3)*(1025 - ctr_3)*(1026 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(1027 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1026 - ctr_3)*(1027 - ctr_3)*(1028 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(1026 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(1026 - ctr_3) + ((1080044550) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((1025 - ctr_3)*(1026 - ctr_3)*(1027 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_11(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 2048; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 2048 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 2048; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(2050 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(2049 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (2051 - ctr_3)*(ctr_2 + 1) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((2050 - ctr_3)*(2051 - ctr_3)*(2052 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (2050 - ctr_3)*(ctr_2 + 1) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (2050 - ctr_3)*(ctr_2 - 1) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (2049 - ctr_3)*(ctr_2 - 1) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(2051 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2050 - ctr_3)*(2051 - ctr_3)*(2052 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(2050 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(2049 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (2051 - ctr_3)*(ctr_2 + 1) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((2050 - ctr_3)*(2051 - ctr_3)*(2052 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (2050 - ctr_3)*(ctr_2 + 1) + ((8615122950) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (2050 - ctr_3)*(ctr_2 - 1) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (2049 - ctr_3)*(ctr_2 - 1) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((2048 - ctr_3)*(2049 - ctr_3)*(2050 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(2051 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2050 - ctr_3)*(2051 - ctr_3)*(2052 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(2050 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(2050 - ctr_3) + ((8615122950) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((2049 - ctr_3)*(2050 - ctr_3)*(2051 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_12(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 4096; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 4096 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4096; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(4098 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(4097 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (4099 - ctr_3)*(ctr_2 + 1) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((4098 - ctr_3)*(4099 - ctr_3)*(4100 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (4098 - ctr_3)*(ctr_2 + 1) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (4098 - ctr_3)*(ctr_2 - 1) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (4097 - ctr_3)*(ctr_2 - 1) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(4099 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4098 - ctr_3)*(4099 - ctr_3)*(4100 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(4098 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(4097 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (4099 - ctr_3)*(ctr_2 + 1) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((4098 - ctr_3)*(4099 - ctr_3)*(4100 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (4098 - ctr_3)*(ctr_2 + 1) + ((68820185094) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (4098 - ctr_3)*(ctr_2 - 1) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (4097 - ctr_3)*(ctr_2 - 1) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((4096 - ctr_3)*(4097 - ctr_3)*(4098 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(4099 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4098 - ctr_3)*(4099 - ctr_3)*(4100 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(4098 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(4098 - ctr_3) + ((68820185094) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((4097 - ctr_3)*(4098 - ctr_3)*(4099 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_13(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 8192; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 8192 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8192; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(8194 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(8193 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (8195 - ctr_3)*(ctr_2 + 1) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((8194 - ctr_3)*(8195 - ctr_3)*(8196 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (8194 - ctr_3)*(ctr_2 + 1) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (8194 - ctr_3)*(ctr_2 - 1) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (8193 - ctr_3)*(ctr_2 - 1) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(8195 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8194 - ctr_3)*(8195 - ctr_3)*(8196 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(8194 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(8193 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (8195 - ctr_3)*(ctr_2 + 1) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((8194 - ctr_3)*(8195 - ctr_3)*(8196 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (8194 - ctr_3)*(ctr_2 + 1) + ((550158557190) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (8194 - ctr_3)*(ctr_2 - 1) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (8193 - ctr_3)*(ctr_2 - 1) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((8192 - ctr_3)*(8193 - ctr_3)*(8194 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(8195 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8194 - ctr_3)*(8195 - ctr_3)*(8196 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(8194 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(8194 - ctr_3) + ((550158557190) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((8193 - ctr_3)*(8194 - ctr_3)*(8195 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_14(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < 16384; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < 16384 - ctr_3; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16384; ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(16386 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(16385 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (16387 - ctr_3)*(ctr_2 + 1) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((16386 - ctr_3)*(16387 - ctr_3)*(16388 - ctr_3)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (16386 - ctr_3)*(ctr_2 + 1) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (16386 - ctr_3)*(ctr_2 - 1) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (16385 - ctr_3)*(ctr_2 - 1) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(16387 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16386 - ctr_3)*(16387 - ctr_3)*(16388 - ctr_3)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(16386 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(16385 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (16387 - ctr_3)*(ctr_2 + 1) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((16386 - ctr_3)*(16387 - ctr_3)*(16388 - ctr_3)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (16386 - ctr_3)*(ctr_2 + 1) + ((4399657304070) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (16386 - ctr_3)*(ctr_2 - 1) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (16385 - ctr_3)*(ctr_2 - 1) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((16384 - ctr_3)*(16385 - ctr_3)*(16386 - ctr_3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(16387 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16386 - ctr_3)*(16387 - ctr_3)*(16388 - ctr_3)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(16386 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(16386 - ctr_3) + ((4399657304070) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((16385 - ctr_3)*(16386 - ctr_3)*(16387 - ctr_3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}

static void apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_any(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_18 = xi_1*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_25 = xi_2*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - 1];
            const double xi_26 = xi_3*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) - 1];
            const double xi_27 = xi_4*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_28 = xi_5*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_29 = xi_6*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_30 = xi_7*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_31 = xi_8*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_32 = xi_9*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_19 = xi_10*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_20 = xi_11*_data_p1CellSrc[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_21 = xi_12*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_22 = xi_13*_data_p1CellSrc[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) + 1];
            const double xi_23 = xi_14*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) + 1];
            const double xi_24 = xi_15*_data_p1CellSrc[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_p1CellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
         }
      }
   }
}


void apply_3D_macrocell_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
    switch( level )
    {
    case 2:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_2(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 3:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_3(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 4:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_4(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 5:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_5(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 6:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_6(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 7:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_7(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 8:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_8(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 9:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_9(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 10:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_10(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 11:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_11(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 12:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_12(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 13:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_13(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    case 14:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_14(_data_p1CellDst, _data_p1CellSrc, p1CellStencil);
        break;
    default:
        apply_3D_macrocell_vertexdof_to_vertexdof_replace_level_any(_data_p1CellDst, _data_p1CellSrc, level, p1CellStencil);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg
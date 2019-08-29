
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "gaussseidel_3D_macrocell_P1_colored_impl_group_1.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void gaussseidel_3D_macrocell_P1_colored_impl_group_1_level_any(double const * RESTRICT const _data_p1CellDst_group_0_const, double * RESTRICT _data_p1CellDst_group_1, double const * RESTRICT const _data_p1CellDst_group_2_const, double const * RESTRICT const _data_p1CellDst_group_3_const, double const * RESTRICT const _data_p1CellDst_group_4_const, double const * RESTRICT const _data_p1CellDst_group_5_const, double const * RESTRICT const _data_p1CellDst_group_6_const, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_1_const, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_115 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_131 = 1 / (xi_115);
   const double xi_116 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_117 = p1CellStencil[{ 1, 0, 0 }];
   const double xi_118 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_119 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_120 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_121 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_122 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_123 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_124 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_125 = p1CellStencil[{ 0, 0, 1 }];
   const double xi_126 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_127 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_128 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_129 = p1CellStencil[{ 0, 1, -1 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level - 1)) - 1; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level - 1)) - 1; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_147 = _data_p1CellRhs_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_133 = -xi_116*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_134 = -xi_117*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_135 = -xi_118*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_136 = -xi_119*_data_p1CellDst_group_2_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_137 = -xi_120*_data_p1CellDst_group_3_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_138 = -xi_121*_data_p1CellDst_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_139 = -xi_122*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_140 = -xi_123*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_141 = -xi_124*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_142 = -xi_125*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_143 = -xi_126*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_144 = -xi_127*_data_p1CellDst_group_6_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_145 = -xi_128*_data_p1CellDst_group_7_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_146 = -xi_129*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_131*(xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146 + xi_147);
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level - 1)) - 1; ctr_1 += 1)
         {
            const double xi_33 = _data_p1CellRhs_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_25 = -xi_116*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_26 = -xi_117*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_27 = -xi_118*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_28 = -xi_119*_data_p1CellDst_group_2_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_29 = -xi_120*_data_p1CellDst_group_3_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_30 = -xi_121*_data_p1CellDst_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_31 = -xi_122*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_32 = -xi_123*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_19 = -xi_124*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_20 = -xi_125*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_21 = -xi_126*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_22 = -xi_127*_data_p1CellDst_group_6_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_23 = -xi_128*_data_p1CellDst_group_7_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_24 = -xi_129*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_131*(xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33);
         }
         // face 3
         for (int ctr_1 = -ctr_2 - ctr_3 + (1 << (level - 1)) - 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level - 1)); ctr_1 += 1)
         {
            const double xi_182 = _data_p1CellRhs_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_168 = -xi_116*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_169 = -xi_117*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_170 = -xi_118*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_171 = -xi_119*_data_p1CellDst_group_2_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_172 = -xi_120*_data_p1CellDst_group_3_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_173 = -xi_121*_data_p1CellDst_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_174 = -xi_122*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_175 = -xi_123*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_176 = -xi_124*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_177 = -xi_125*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_178 = -xi_126*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_179 = -xi_127*_data_p1CellDst_group_6_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_180 = -xi_128*_data_p1CellDst_group_7_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_181 = -xi_129*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_131*(xi_168 + xi_169 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182);
         }
      }
      for (int ctr_2 = -ctr_3 + (1 << (level - 1)) - 1; ctr_2 < -ctr_3 + (1 << (level - 1)); ctr_2 += 1)
      {
         // edge 5
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_104 = _data_p1CellRhs_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_90 = -xi_116*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_91 = -xi_117*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_92 = -xi_118*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_93 = -xi_119*_data_p1CellDst_group_2_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_94 = -xi_120*_data_p1CellDst_group_3_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_95 = -xi_121*_data_p1CellDst_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_96 = -xi_122*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_97 = -xi_123*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_98 = -xi_124*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_99 = -xi_125*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_100 = -xi_126*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_101 = -xi_127*_data_p1CellDst_group_6_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_102 = -xi_128*_data_p1CellDst_group_7_const[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_103 = -xi_129*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_131*(xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99);
         }
      }
   }
}


void gaussseidel_3D_macrocell_P1_colored_impl_group_1(double const * RESTRICT const _data_p1CellDst_group_0_const, double * RESTRICT _data_p1CellDst_group_1, double const * RESTRICT const _data_p1CellDst_group_2_const, double const * RESTRICT const _data_p1CellDst_group_3_const, double const * RESTRICT const _data_p1CellDst_group_4_const, double const * RESTRICT const _data_p1CellDst_group_5_const, double const * RESTRICT const _data_p1CellDst_group_6_const, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_1_const, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
    switch( level )
    {

    default:
        gaussseidel_3D_macrocell_P1_colored_impl_group_1_level_any(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_1_const, level, p1CellStencil);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg
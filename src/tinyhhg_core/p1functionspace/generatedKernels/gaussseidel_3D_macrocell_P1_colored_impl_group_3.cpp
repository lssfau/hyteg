
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "gaussseidel_3D_macrocell_P1_colored_impl_group_3.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void gaussseidel_3D_macrocell_P1_colored_impl_group_3_level_any(double const * RESTRICT const _data_p1CellDst_group_0_const, double const * RESTRICT const _data_p1CellDst_group_1_const, double const * RESTRICT const _data_p1CellDst_group_2_const, double * RESTRICT _data_p1CellDst_group_3, double const * RESTRICT const _data_p1CellDst_group_4_const, double const * RESTRICT const _data_p1CellDst_group_5_const, double const * RESTRICT const _data_p1CellDst_group_6_const, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_3_const, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
   const double xi_64 = p1CellStencil[{ 0, 0, 0 }];
   const double xi_80 = 1 / (xi_64);
   const double xi_65 = p1CellStencil[{ -1, 1, 0 }];
   const double xi_66 = p1CellStencil[{ 1, -1, 0 }];
   const double xi_67 = p1CellStencil[{ 0, -1, 0 }];
   const double xi_68 = p1CellStencil[{ 0, 1, 0 }];
   const double xi_69 = p1CellStencil[{ -1, 0, 0 }];
   const double xi_70 = p1CellStencil[{ 1, 0, 0 }];
   const double xi_71 = p1CellStencil[{ -1, 1, -1 }];
   const double xi_72 = p1CellStencil[{ 1, -1, 1 }];
   const double xi_73 = p1CellStencil[{ 0, -1, 1 }];
   const double xi_74 = p1CellStencil[{ 0, 1, -1 }];
   const double xi_75 = p1CellStencil[{ -1, 0, 1 }];
   const double xi_76 = p1CellStencil[{ 1, 0, -1 }];
   const double xi_77 = p1CellStencil[{ 0, 0, -1 }];
   const double xi_78 = p1CellStencil[{ 0, 0, 1 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level - 1)) - 1; ctr_3 += 1)
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // edge 3
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_96 = _data_p1CellRhs_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_82 = -xi_65*_data_p1CellDst_group_0_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_83 = -xi_66*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_84 = -xi_67*_data_p1CellDst_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_85 = -xi_68*_data_p1CellDst_group_1_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_86 = -xi_69*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_87 = -xi_70*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_88 = -xi_71*_data_p1CellDst_group_4_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_89 = -xi_72*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_90 = -xi_73*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_91 = -xi_74*_data_p1CellDst_group_5_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_92 = -xi_75*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_93 = -xi_76*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_94 = -xi_77*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_95 = -xi_78*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_80*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96);
         }
         // face 1
         for (int ctr_1 = 1; ctr_1 < -ctr_3 + (1 << (level - 1)) - 1; ctr_1 += 1)
         {
            const double xi_143 = _data_p1CellRhs_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_129 = -xi_65*_data_p1CellDst_group_0_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_130 = -xi_66*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_131 = -xi_67*_data_p1CellDst_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_132 = -xi_68*_data_p1CellDst_group_1_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_133 = -xi_69*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_134 = -xi_70*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_135 = -xi_71*_data_p1CellDst_group_4_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_136 = -xi_72*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_137 = -xi_73*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_138 = -xi_74*_data_p1CellDst_group_5_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_139 = -xi_75*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_140 = -xi_76*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_141 = -xi_77*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_142 = -xi_78*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_80*(xi_129 + xi_130 + xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143);
         }
      }
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level - 1)) - 1; ctr_2 += 1)
      {
         // face 2
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_178 = _data_p1CellRhs_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_164 = -xi_65*_data_p1CellDst_group_0_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_165 = -xi_66*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_166 = -xi_67*_data_p1CellDst_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_167 = -xi_68*_data_p1CellDst_group_1_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_168 = -xi_69*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_169 = -xi_70*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_170 = -xi_71*_data_p1CellDst_group_4_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_171 = -xi_72*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_172 = -xi_73*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_173 = -xi_74*_data_p1CellDst_group_5_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_174 = -xi_75*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_175 = -xi_76*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_176 = -xi_77*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_177 = -xi_78*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_80*(xi_164 + xi_165 + xi_166 + xi_167 + xi_168 + xi_169 + xi_170 + xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178);
         }
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level - 1)) - 1; ctr_1 += 1)
         {
            const double xi_33 = _data_p1CellRhs_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_25 = -xi_65*_data_p1CellDst_group_0_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_26 = -xi_66*_data_p1CellDst_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_27 = -xi_67*_data_p1CellDst_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_28 = -xi_68*_data_p1CellDst_group_1_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_29 = -xi_69*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_30 = -xi_70*_data_p1CellDst_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_31 = -xi_71*_data_p1CellDst_group_4_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_32 = -xi_72*_data_p1CellDst_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) + 1];
            const double xi_19 = -xi_73*_data_p1CellDst_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_20 = -xi_74*_data_p1CellDst_group_5_const[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level - 1)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
            const double xi_21 = -xi_75*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
            const double xi_22 = -xi_76*_data_p1CellDst_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6)) + 1];
            const double xi_23 = -xi_77*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            const double xi_24 = -xi_78*_data_p1CellDst_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
            _data_p1CellDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_80*(xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33);
         }
      }
   }
}


void gaussseidel_3D_macrocell_P1_colored_impl_group_3(double const * RESTRICT const _data_p1CellDst_group_0_const, double const * RESTRICT const _data_p1CellDst_group_1_const, double const * RESTRICT const _data_p1CellDst_group_2_const, double * RESTRICT _data_p1CellDst_group_3, double const * RESTRICT const _data_p1CellDst_group_4_const, double const * RESTRICT const _data_p1CellDst_group_5_const, double const * RESTRICT const _data_p1CellDst_group_6_const, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_3_const, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil)
{
    switch( level )
    {

    default:
        gaussseidel_3D_macrocell_P1_colored_impl_group_3_level_any(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_3_const, level, p1CellStencil);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macrocell_P1_colored.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P1_colored(double * RESTRICT _data_p1CellDst_group_0, double const * RESTRICT const _data_p1CellDst_group_0_const, double * RESTRICT _data_p1CellDst_group_1, double const * RESTRICT const _data_p1CellDst_group_1_const, double * RESTRICT _data_p1CellDst_group_2, double const * RESTRICT const _data_p1CellDst_group_2_const, double * RESTRICT _data_p1CellDst_group_3, double const * RESTRICT const _data_p1CellDst_group_3_const, double * RESTRICT _data_p1CellDst_group_4, double const * RESTRICT const _data_p1CellDst_group_4_const, double * RESTRICT _data_p1CellDst_group_5, double const * RESTRICT const _data_p1CellDst_group_5_const, double * RESTRICT _data_p1CellDst_group_6, double const * RESTRICT const _data_p1CellDst_group_6_const, double * RESTRICT _data_p1CellDst_group_7, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_0_const, double const * RESTRICT const _data_p1CellRhs_group_1_const, double const * RESTRICT const _data_p1CellRhs_group_2_const, double const * RESTRICT const _data_p1CellRhs_group_3_const, double const * RESTRICT const _data_p1CellRhs_group_4_const, double const * RESTRICT const _data_p1CellRhs_group_5_const, double const * RESTRICT const _data_p1CellRhs_group_6_const, double const * RESTRICT const _data_p1CellRhs_group_7_const, int64_t color_group, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil, double relax)
{
   if (((color_group) == (0)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_0(_data_p1CellDst_group_0, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_0_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (1)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_1(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_1_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (2)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_2(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_2_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (3)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_3(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_3_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (4)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_4(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_4_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (5)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_5(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_5_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (6)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_6(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6, _data_p1CellDst_group_7_const, _data_p1CellRhs_group_6_const, level, p1CellStencil, relax);
      
      return;
   } 
   if (((color_group) == (7)))
   {
      
      sor_3D_macrocell_P1_colored_impl_group_7(_data_p1CellDst_group_0_const, _data_p1CellDst_group_1_const, _data_p1CellDst_group_2_const, _data_p1CellDst_group_3_const, _data_p1CellDst_group_4_const, _data_p1CellDst_group_5_const, _data_p1CellDst_group_6_const, _data_p1CellDst_group_7, _data_p1CellRhs_group_7_const, level, p1CellStencil, relax);
      
      return;
   } 
}


} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P1_colored(double * RESTRICT _data_p1CellDst_group_0, double const * RESTRICT const _data_p1CellDst_group_0_const, double * RESTRICT _data_p1CellDst_group_1, double const * RESTRICT const _data_p1CellDst_group_1_const, double * RESTRICT _data_p1CellDst_group_2, double const * RESTRICT const _data_p1CellDst_group_2_const, double * RESTRICT _data_p1CellDst_group_3, double const * RESTRICT const _data_p1CellDst_group_3_const, double * RESTRICT _data_p1CellDst_group_4, double const * RESTRICT const _data_p1CellDst_group_4_const, double * RESTRICT _data_p1CellDst_group_5, double const * RESTRICT const _data_p1CellDst_group_5_const, double * RESTRICT _data_p1CellDst_group_6, double const * RESTRICT const _data_p1CellDst_group_6_const, double * RESTRICT _data_p1CellDst_group_7, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_0_const, double const * RESTRICT const _data_p1CellRhs_group_1_const, double const * RESTRICT const _data_p1CellRhs_group_2_const, double const * RESTRICT const _data_p1CellRhs_group_3_const, double const * RESTRICT const _data_p1CellRhs_group_4_const, double const * RESTRICT const _data_p1CellRhs_group_5_const, double const * RESTRICT const _data_p1CellRhs_group_6_const, double const * RESTRICT const _data_p1CellRhs_group_7_const, int64_t color_group, int32_t level, std::map< hhg::indexing::IndexIncrement, double > p1CellStencil, double relax);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg
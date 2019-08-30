
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P1(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellRhs, int32_t level, std::map< hyteg::indexing::IndexIncrement, double > p1CellStencil, double relax);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg
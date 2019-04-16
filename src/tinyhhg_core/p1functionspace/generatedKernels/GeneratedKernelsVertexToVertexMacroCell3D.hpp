
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "core/Macros.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace vertexdof {
namespace macrocell {
namespace generated {

void apply_3D_macrocell_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, double const * const _data_p1CellStencil, int64_t level);

void apply_3D_macrocell_vertexdof_to_vertexdof_add(double * RESTRICT _data_p1CellDstAdd, double const * RESTRICT const _data_p1CellSrcAdd, double const * const _data_p1CellStencil, int64_t level);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg
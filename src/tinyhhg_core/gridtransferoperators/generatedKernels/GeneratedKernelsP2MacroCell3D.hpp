
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
namespace P2 {
namespace macrocell {
namespace generated {

void restrict_3D_macrocell_P2_update_vertexdofs(double const * RESTRICT const _data_edgeFineSrc_X, double const * RESTRICT const _data_edgeFineSrc_XY, double const * RESTRICT const _data_edgeFineSrc_XYZ, double const * RESTRICT const _data_edgeFineSrc_XZ, double const * RESTRICT const _data_edgeFineSrc_Y, double const * RESTRICT const _data_edgeFineSrc_YZ, double const * RESTRICT const _data_edgeFineSrc_Z, double * RESTRICT _data_vertexCoarseDst, double const * RESTRICT const _data_vertexFineSrc, int64_t coarse_level);

void restrict_3D_macrocell_P2_update_edgedofs(double * RESTRICT _data_edgeCoarseDst_X, double * RESTRICT _data_edgeCoarseDst_XY, double * RESTRICT _data_edgeCoarseDst_XYZ, double * RESTRICT _data_edgeCoarseDst_XZ, double * RESTRICT _data_edgeCoarseDst_Y, double * RESTRICT _data_edgeCoarseDst_YZ, double * RESTRICT _data_edgeCoarseDst_Z, double const * RESTRICT const _data_edgeFineSrc_X, double const * RESTRICT const _data_edgeFineSrc_XY, double const * RESTRICT const _data_edgeFineSrc_XYZ, double const * RESTRICT const _data_edgeFineSrc_XZ, double const * RESTRICT const _data_edgeFineSrc_Y, double const * RESTRICT const _data_edgeFineSrc_YZ, double const * RESTRICT const _data_edgeFineSrc_Z, double const * RESTRICT const _data_vertexFineSrc, int64_t coarse_level);

void prolongate_3D_macrocell_P2_push_from_vertexdofs(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_XYZ, double * RESTRICT _data_edgeFineDst_XZ, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_edgeFineDst_YZ, double * RESTRICT _data_edgeFineDst_Z, double const * RESTRICT const _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int64_t coarse_level);

void prolongate_3D_macrocell_P2_push_from_edgedofs(double const * RESTRICT const _data_edgeCoarseSrc_X, double const * RESTRICT const _data_edgeCoarseSrc_XY, double const * RESTRICT const _data_edgeCoarseSrc_XYZ, double const * RESTRICT const _data_edgeCoarseSrc_XZ, double const * RESTRICT const _data_edgeCoarseSrc_Y, double const * RESTRICT const _data_edgeCoarseSrc_YZ, double const * RESTRICT const _data_edgeCoarseSrc_Z, double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_XYZ, double * RESTRICT _data_edgeFineDst_XZ, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_edgeFineDst_YZ, double * RESTRICT _data_edgeFineDst_Z, double * RESTRICT _data_vertexFineDst, int64_t coarse_level);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hhg
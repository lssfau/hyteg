#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

void applyReplace( double*          fd_edgeFaceDst,
                   double*          fd_edgeFaceSrc,
                   double*          fd_edgeFaceStencil,
                   walberla::uint_t level );

void applyAdd( double* fd_edgeFaceDst, double* fd_edgeFaceSrc, double* fd_edgeFaceStencil, walberla::uint_t level );

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
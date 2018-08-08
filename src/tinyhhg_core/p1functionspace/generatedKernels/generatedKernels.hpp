#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

void applyReplace( double*                fd_p1FaceDst,
                   const double*          fd_p1FaceSrc,
                   const double*          fd_p1FaceStencil,
                   const walberla::uint_t level );

void applyAdd( double* fd_p1FaceDst, const double* fd_p1FaceSrc, const double* fd_p1FaceStencil, const walberla::uint_t level );

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
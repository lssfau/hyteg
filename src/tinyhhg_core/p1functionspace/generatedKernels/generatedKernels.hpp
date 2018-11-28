#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

void applyReplace( double*                fd_p1FaceDst,
                   double*          fd_p1FaceSrc,
                   double*          fd_p1FaceStencil,
                   walberla::uint_t level );

void applyAdd( double* fd_p1FaceDst, double* fd_p1FaceSrc, double* fd_p1FaceStencil, walberla::uint_t level );

void gaussSeidel( double * fd_p1FaceDst, double * fd_p1FaceRhs, double * fd_p1FaceStencil, walberla::uint_t level );

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
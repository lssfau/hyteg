#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void applyFaceReplace( double*                fd_p1FaceDst,
                       double*          fd_edgeFaceSrc,
                       double*          fd_edgeFaceStencil,
                       walberla::uint_t level );

void applyFaceAdd( double* fd_p1FaceDst, double* fd_edgeFaceSrc, double* fd_edgeFaceStencil, walberla::uint_t level );

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
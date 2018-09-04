#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

void applyFaceReplace( double*          fd_edgeFaceDst,
                       double*          fd_vertexFaceSrc,
                       double*          fd_vertexToEdgeFaceStencil,
                       walberla::uint_t level );

void applyFaceAdd( double*          fd_edgeFaceDst,
                   double*          fd_vertexFaceSrc,
                   double*          fd_vertexToEdgeFaceStencil,
                   walberla::uint_t level );


} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
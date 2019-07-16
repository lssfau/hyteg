#pragma once

#include "core/DataTypes.h"
#include "core/Macros.h"

void apply_2d_vertex_to_vertex( double* WALBERLA_RESTRICT dstPtr,
                                double const* WALBERLA_RESTRICT const srcPtr,
                                double const* WALBERLA_RESTRICT const stencilPtr,
                                const walberla::uint_t                          rowsize );
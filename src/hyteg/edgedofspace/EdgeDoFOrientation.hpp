#pragma once

#include "core/DataTypes.h"

namespace hyteg {
namespace edgedof {

/// Mapping of X,Y,Z coordinates to uint_t
/// located in seperate file to reduce dependencies in generated kernels
enum class EdgeDoFOrientation : walberla::uint_t
{
   X,
   Y,
   Z,
   XY,
   XZ,
   YZ,
   XYZ,
   INVALID,
};
}// namespace edgedof
}// namespace hyteg
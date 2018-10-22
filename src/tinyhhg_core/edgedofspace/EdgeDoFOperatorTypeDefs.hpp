#pragma once

#include <map>

#include "core/DataTypes.h"

namespace hhg {

namespace indexing {
class IndexIncrement;
}

namespace edgedof {

using walberla::real_t;
using walberla::uint_t;

enum class EdgeDoFOrientation : walberla::uint_t;

namespace macroedge {

/// map[neighborCellID][centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< uint_t,
                  std::map< edgedof::EdgeDoFOrientation,
                            std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > >
    StencilMap_T;

} // namespace macroedge

namespace macroface {

/// map[neighborCellID][centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< uint_t,
                  std::map< edgedof::EdgeDoFOrientation,
                            std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > >
    StencilMap_T;

} // namespace macroface

namespace macrocell {

/// map[centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< edgedof::EdgeDoFOrientation,
                  std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > >
    StencilMap_T;

} // namespace macrocell

} // namespace edgedof
} // namespace hhg
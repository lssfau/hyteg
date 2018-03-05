#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"

namespace hhg {
namespace DGEdge{

using walberla::uint_t;

constexpr inline uint_t indexDGFaceFromVertex( const uint_t & level, const uint_t pos, stencilDirection dir ) {
  return hhg::BubbleEdge::indexFaceFromVertex( level, pos, dir );
}




}// namesapce DGEdge
}// namesspace hhg
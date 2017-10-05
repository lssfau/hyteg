#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"

namespace hhg {
namespace DGEdge{

using walberla::uint_t;

template<uint_t Level>
constexpr inline uint_t indexDGFaceFromVertex(const uint_t pos, stencilDirection dir){
  return 1;
  //hhg::BubbleEdge::EdgeCoordsVertex::index<3>(pos,);
}




}// namesapce DGEdge
}// namesspace hhg
#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

namespace hhg {
namespace DgFace{

using walberla::uint_t;

template<size_t Level>
constexpr inline uint_t indexDGcellFromVertex( const uint_t col, const uint_t row, const stencilDirection dir){
  return BubbleFace::indexFaceFromVertex<Level>(col, row, dir);
}

template<size_t Level>
constexpr inline uint_t indexDGcellFromGrayDGCell( const uint_t col, const uint_t row, const stencilDirection dir){
  return BubbleFace::indexFaceFromVertex<Level>(col - 1, row, dir);
}




}// namespace DgFace

}// namespace hhg
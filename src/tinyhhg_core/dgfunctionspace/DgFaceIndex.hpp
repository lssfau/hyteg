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
  typedef hhg::stencilDirection sD;
  switch(dir){
    case sD::CELL_BLUE_S:
      return BubbleFace::indexFaceFromVertex<Level>(col, row, sD::CELL_BLUE_SE);
    case sD::CELL_BLUE_E:
      return BubbleFace::indexFaceFromVertex<Level>(col + 1, row, sD::CELL_BLUE_NW);
    case sD::CELL_BLUE_W:
      return BubbleFace::indexFaceFromVertex<Level>(col, row, sD::CELL_BLUE_NW);
  }
  WALBERLA_ASSERT(false);
  return std::numeric_limits<size_t>::max();
}

template<size_t Level>
constexpr inline uint_t indexDGcellFromBlueDGCell( const uint_t col, const uint_t row, const stencilDirection dir){
  typedef hhg::stencilDirection sD;
  switch(dir){
    case sD::CELL_GRAY_E:
      return BubbleFace::indexFaceFromVertex<Level>(col + 1, row, sD::CELL_GRAY_NE);
    case sD::CELL_GRAY_N:
      return BubbleFace::indexFaceFromVertex<Level>(col, row + 1, sD::CELL_GRAY_NE);
    case sD::CELL_GRAY_W:
      return BubbleFace::indexFaceFromVertex<Level>(col, row, sD::CELL_GRAY_NE);
  }
  WALBERLA_ASSERT(false);
  return std::numeric_limits<size_t>::max();
}



}// namespace DgFace

}// namespace hhg
#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

namespace hhg {
namespace DGFace{

using walberla::uint_t;

template<size_t Level>
constexpr inline uint_t indexDGFaceFromVertex(const uint_t col, const uint_t row, const stencilDirection dir){
  return BubbleFace::indexFaceFromVertex( Level, col, row, dir );
}

constexpr std::array<hhg::stencilDirection ,3> grayDGfaceneighbors =
  {{stencilDirection::CELL_BLUE_S, stencilDirection::CELL_BLUE_E, stencilDirection::CELL_BLUE_W}};

template<size_t Level>
constexpr inline uint_t indexDGFaceFromGrayDGface(const uint_t col, const uint_t row, const stencilDirection dir){
  typedef hhg::stencilDirection sD;
  switch(dir){
    case sD::CELL_BLUE_S:
      return BubbleFace::indexFaceFromVertex( Level, col, row, sD::CELL_BLUE_SE );
    case sD::CELL_BLUE_E:
      return BubbleFace::indexFaceFromVertex( Level, col + 1, row, sD::CELL_BLUE_NW );
    case sD::CELL_BLUE_W:
      return BubbleFace::indexFaceFromVertex( Level, col, row, sD::CELL_BLUE_NW );
  }
  WALBERLA_ASSERT(false);
  return std::numeric_limits<size_t>::max();
}

constexpr std::array<hhg::stencilDirection ,3> blueDGfaceneighbors =
  {{stencilDirection::CELL_GRAY_E, stencilDirection::CELL_GRAY_N, stencilDirection::CELL_GRAY_W}};

template<size_t Level>
constexpr inline uint_t indexDGFaceFromBlueDGface(const uint_t col, const uint_t row, const stencilDirection dir){
  typedef hhg::stencilDirection sD;
  switch(dir){
    case sD::CELL_GRAY_E:
      return BubbleFace::indexFaceFromVertex( Level, col + 1, row, sD::CELL_GRAY_NE );
    case sD::CELL_GRAY_N:
      return BubbleFace::indexFaceFromVertex( Level, col, row + 1, sD::CELL_GRAY_NE );
    case sD::CELL_GRAY_W:
      return BubbleFace::indexFaceFromVertex( Level, col, row, sD::CELL_GRAY_NE );
  }
  WALBERLA_ASSERT(false);
  return std::numeric_limits<size_t>::max();
}

//we can use the same indexIterator as in Bubble
using BubbleFace::indexIterator;

}// namespace DGFace

}// namespace hhg
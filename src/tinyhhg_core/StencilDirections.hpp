#pragma once

#include "core/DataTypes.h"

using walberla::uint_t;

namespace hhg {

enum class stencilDirection : uint_t {
  VERTEX_C = 0,
  VERTEX_S = 1,
  VERTEX_SE = 2,
  VERTEX_E = 3,
  VERTEX_NE = 4,
  VERTEX_N = 5,
  VERTEX_NW = 6,
  VERTEX_W = 7,
  VERTEX_SW = 8,

  CELL_BLUE_C = 9,
  CELL_BLUE_S = 10,
  CELL_BLUE_SE = 11,
  CELL_BLUE_E = 12,
  CELL_BLUE_NE = 13,
  CELL_BLUE_N = 14,
  CELL_BLUE_NW = 15,
  CELL_BLUE_W = 16,
  CELL_BLUE_SW = 17,

  CELL_GRAY_C = 18,
  CELL_GRAY_S = 19,
  CELL_GRAY_SE = 20,
  CELL_GRAY_E = 21,
  CELL_GRAY_NE = 24,
  CELL_GRAY_N = 25,
  CELL_GRAY_NW = 26,
  CELL_GRAY_W = 27,
  CELL_GRAY_SW = 28,

  EDGE_HO_C,
  EDGE_HO_W,
  EDGE_HO_E,
  EDGE_HO_NW,
  EDGE_HO_SE,

  EDGE_VE_N,
  EDGE_VE_S,
  EDGE_VE_NW,
  EDGE_VE_SE,

  EDGE_DI_N,
  EDGE_DI_S,
  EDGE_DI_NW,
  EDGE_DI_NE,
  EDGE_DI_SW,
  EDGE_DI_SE,

};

}//namespace hhg

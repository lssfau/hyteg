#pragma once

#include "core/DataTypes.h"

using walberla::uint_t;

namespace hhg {

enum class stencilDirection : uint_t {
  VERTEX_C = 0,
  VERTEX_S,
  VERTEX_SE,
  VERTEX_E,
  VERTEX_NE,
  VERTEX_N,
  VERTEX_NW,
  VERTEX_W,
  VERTEX_SW,

  VERTEX_BC,
  VERTEX_BW,
  VERTEX_BS,
  VERTEX_BSW,

  VERTEX_FC,
  VERTEX_FE,
  VERTEX_FN,
  VERTEX_FNE,

  CELL_BLUE_C,
  CELL_BLUE_S,
  CELL_BLUE_SE,
  CELL_BLUE_E,
  CELL_BLUE_NE,
  CELL_BLUE_N,
  CELL_BLUE_NW,
  CELL_BLUE_W,
  CELL_BLUE_SW,

  CELL_GRAY_C,
  CELL_GRAY_S,
  CELL_GRAY_SE,
  CELL_GRAY_E,
  CELL_GRAY_NE,
  CELL_GRAY_N,
  CELL_GRAY_NW,
  CELL_GRAY_W,
  CELL_GRAY_SW,

  EDGE_HO_C,
  EDGE_HO_N,
  EDGE_HO_S,
  EDGE_HO_W,
  EDGE_HO_E,
  EDGE_HO_NW,
  EDGE_HO_SE,

  EDGE_VE_C,
  EDGE_VE_N,
  EDGE_VE_S,
  EDGE_VE_W,
  EDGE_VE_E,
  EDGE_VE_NW,
  EDGE_VE_SE,

  EDGE_DI_C,
  EDGE_DI_N,
  EDGE_DI_S,
  EDGE_DI_W,
  EDGE_DI_E,
  EDGE_DI_NW,
  EDGE_DI_NE,
  EDGE_DI_SW,
  EDGE_DI_SE,

};

inline bool isHorizontalEdge( const stencilDirection & dir )
{
  switch ( dir )
  {
  case stencilDirection::EDGE_HO_C:
  case stencilDirection::EDGE_HO_N:
  case stencilDirection::EDGE_HO_S:
  case stencilDirection::EDGE_HO_W:
  case stencilDirection::EDGE_HO_E:
  case stencilDirection::EDGE_HO_NW:
  case stencilDirection::EDGE_HO_SE:
    return true;
    break;
  default:
    return false;
    break;
  }
}

inline bool isDiagonalEdge( const stencilDirection & dir )
{
  switch ( dir )
  {
  case stencilDirection::EDGE_DI_C:
  case stencilDirection::EDGE_DI_N:
  case stencilDirection::EDGE_DI_S:
  case stencilDirection::EDGE_DI_W:
  case stencilDirection::EDGE_DI_E:
  case stencilDirection::EDGE_DI_NW:
  case stencilDirection::EDGE_DI_NE:
  case stencilDirection::EDGE_DI_SW:
  case stencilDirection::EDGE_DI_SE:
    return true;
    break;
  default:
    return false;
    break;
  }
}

inline bool isVerticalEdge( const stencilDirection & dir )
{
  switch ( dir )
  {
  case stencilDirection::EDGE_VE_C:
  case stencilDirection::EDGE_VE_N:
  case stencilDirection::EDGE_VE_S:
  case stencilDirection::EDGE_VE_W:
  case stencilDirection::EDGE_VE_E:
  case stencilDirection::EDGE_VE_NW:
  case stencilDirection::EDGE_VE_SE:
    return true;
    break;
  default:
    return false;
    break;
  }
}

}//namespace hhg

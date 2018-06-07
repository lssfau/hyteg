#pragma once

#include <map>
#include "core/DataTypes.h"

#include <map>

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

  VERTEX_TC,
  VERTEX_TS,
  VERTEX_TSE,
  VERTEX_TE,
  VERTEX_TNE,
  VERTEX_TN,
  VERTEX_TNW,
  VERTEX_TW,
  VERTEX_TSW,

  VERTEX_BC,
  VERTEX_BS,
  VERTEX_BSE,
  VERTEX_BE,
  VERTEX_BNE,
  VERTEX_BN,
  VERTEX_BNW,
  VERTEX_BW,
  VERTEX_BSW,

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


static std::map<stencilDirection,std::string> stencilDirectionToStr {

  std::make_pair( stencilDirection::VERTEX_C,    "VERTEX_C" ),
  std::make_pair( stencilDirection::VERTEX_S,    "VERTEX_S" ),
  std::make_pair( stencilDirection::VERTEX_SE,   "VERTEX_SE" ),
  std::make_pair( stencilDirection::VERTEX_E,    "VERTEX_E" ),
  std::make_pair( stencilDirection::VERTEX_NE,   "VERTEX_NE" ),
  std::make_pair( stencilDirection::VERTEX_N,    "VERTEX_N" ),
  std::make_pair( stencilDirection::VERTEX_NW,   "VERTEX_NW" ),
  std::make_pair( stencilDirection::VERTEX_W,    "VERTEX_W" ),
  std::make_pair( stencilDirection::VERTEX_SW,   "VERTEX_SW" ),

  std::make_pair( stencilDirection::VERTEX_TC,    "VERTEX_TC" ),
  std::make_pair( stencilDirection::VERTEX_TS,    "VERTEX_TS" ),
  std::make_pair( stencilDirection::VERTEX_TSE,   "VERTEX_TSE" ),
  std::make_pair( stencilDirection::VERTEX_TE,    "VERTEX_TE" ),
  std::make_pair( stencilDirection::VERTEX_TNE,   "VERTEX_TNE" ),
  std::make_pair( stencilDirection::VERTEX_TN,    "VERTEX_TN" ),
  std::make_pair( stencilDirection::VERTEX_TNW,   "VERTEX_TNW" ),
  std::make_pair( stencilDirection::VERTEX_TW,    "VERTEX_TW" ),
  std::make_pair( stencilDirection::VERTEX_TSW,   "VERTEX_TSW" ),

  std::make_pair( stencilDirection::VERTEX_BC,    "VERTEX_BC" ),
  std::make_pair( stencilDirection::VERTEX_BS,    "VERTEX_BS" ),
  std::make_pair( stencilDirection::VERTEX_BSE,   "VERTEX_BSE" ),
  std::make_pair( stencilDirection::VERTEX_BE,    "VERTEX_BE" ),
  std::make_pair( stencilDirection::VERTEX_BNE,   "VERTEX_BNE" ),
  std::make_pair( stencilDirection::VERTEX_BN,    "VERTEX_BN" ),
  std::make_pair( stencilDirection::VERTEX_BNW,   "VERTEX_BNW" ),
  std::make_pair( stencilDirection::VERTEX_BW,    "VERTEX_BW" ),
  std::make_pair( stencilDirection::VERTEX_BSW,   "VERTEX_BSW" ),

  std::make_pair( stencilDirection::CELL_BLUE_C, "CELL_BLUE_C" ),
  std::make_pair( stencilDirection::CELL_BLUE_S, "CELL_BLUE_S" ),
  std::make_pair( stencilDirection::CELL_BLUE_SE,"CELL_BLUE_SE" ),
  std::make_pair( stencilDirection::CELL_BLUE_E, "CELL_BLUE_E" ),
  std::make_pair( stencilDirection::CELL_BLUE_NE,"CELL_BLUE_NE" ),
  std::make_pair( stencilDirection::CELL_BLUE_N, "CELL_BLUE_N" ),
  std::make_pair( stencilDirection::CELL_BLUE_NW,"CELL_BLUE_NW" ),
  std::make_pair( stencilDirection::CELL_BLUE_W, "CELL_BLUE_W" ),
  std::make_pair( stencilDirection::CELL_BLUE_SW,"CELL_BLUE_SW" ),

  std::make_pair( stencilDirection::CELL_GRAY_C, "CELL_GRAY_C" ),
  std::make_pair( stencilDirection::CELL_GRAY_S, "CELL_GRAY_S" ),
  std::make_pair( stencilDirection::CELL_GRAY_SE,"CELL_GRAY_SE" ),
  std::make_pair( stencilDirection::CELL_GRAY_E, "CELL_GRAY_E" ),
  std::make_pair( stencilDirection::CELL_GRAY_NE,"CELL_GRAY_NE" ),
  std::make_pair( stencilDirection::CELL_GRAY_N, "CELL_GRAY_N" ),
  std::make_pair( stencilDirection::CELL_GRAY_NW,"CELL_GRAY_NW" ),
  std::make_pair( stencilDirection::CELL_GRAY_W, "CELL_GRAY_W" ),
  std::make_pair( stencilDirection::CELL_GRAY_SW,"CELL_GRAY_SW" ),

  std::make_pair( stencilDirection::EDGE_HO_C,   "EDGE_HO_C" ),
  std::make_pair( stencilDirection::EDGE_HO_N,   "EDGE_HO_N" ),
  std::make_pair( stencilDirection::EDGE_HO_S,   "EDGE_HO_S" ),
  std::make_pair( stencilDirection::EDGE_HO_W,   "EDGE_HO_W" ),
  std::make_pair( stencilDirection::EDGE_HO_E,   "EDGE_HO_E" ),
  std::make_pair( stencilDirection::EDGE_HO_NW,  "EDGE_HO_NW" ),
  std::make_pair( stencilDirection::EDGE_HO_SE,  "EDGE_HO_SE" ),

  std::make_pair( stencilDirection::EDGE_VE_C,   "EDGE_VE_C" ),
  std::make_pair( stencilDirection::EDGE_VE_N,   "EDGE_VE_N" ),
  std::make_pair( stencilDirection::EDGE_VE_S,   "EDGE_VE_S" ),
  std::make_pair( stencilDirection::EDGE_VE_W,   "EDGE_VE_W" ),
  std::make_pair( stencilDirection::EDGE_VE_E,   "EDGE_VE_E" ),
  std::make_pair( stencilDirection::EDGE_VE_NW,  "EDGE_VE_NW" ),
  std::make_pair( stencilDirection::EDGE_VE_SE,  "EDGE_VE_SE" ),

  std::make_pair( stencilDirection::EDGE_DI_C,   "EDGE_DI_C" ),
  std::make_pair( stencilDirection::EDGE_DI_N,   "EDGE_DI_N" ),
  std::make_pair( stencilDirection::EDGE_DI_S,   "EDGE_DI_S" ),
  std::make_pair( stencilDirection::EDGE_DI_W,   "EDGE_DI_W" ),
  std::make_pair( stencilDirection::EDGE_DI_E,   "EDGE_DI_E" ),
  std::make_pair( stencilDirection::EDGE_DI_NW,  "EDGE_DI_NW" ),
  std::make_pair( stencilDirection::EDGE_DI_NE,  "EDGE_DI_NE" ),
  std::make_pair( stencilDirection::EDGE_DI_SW,  "EDGE_DI_SW" ),
  std::make_pair( stencilDirection::EDGE_DI_SE,  "EDGE_DI_SE" )

};


inline bool isVertex( const stencilDirection & dir )
{
  switch ( dir )
  {
    case stencilDirection::VERTEX_C:
    case stencilDirection::VERTEX_S:
    case stencilDirection::VERTEX_SE:
    case stencilDirection::VERTEX_E:
    case stencilDirection::VERTEX_NE:
    case stencilDirection::VERTEX_N:
    case stencilDirection::VERTEX_NW:
    case stencilDirection::VERTEX_W:
    case stencilDirection::VERTEX_SW:
    
    case stencilDirection::VERTEX_TC:
    case stencilDirection::VERTEX_TS:
    case stencilDirection::VERTEX_TSE:
    case stencilDirection::VERTEX_TE:
    case stencilDirection::VERTEX_TNE:
    case stencilDirection::VERTEX_TN:
    case stencilDirection::VERTEX_TNW:
    case stencilDirection::VERTEX_TW:
    case stencilDirection::VERTEX_TSW:

    case stencilDirection::VERTEX_BC:
    case stencilDirection::VERTEX_BS:
    case stencilDirection::VERTEX_BSE:
    case stencilDirection::VERTEX_BE:
    case stencilDirection::VERTEX_BNE:
    case stencilDirection::VERTEX_BN:
    case stencilDirection::VERTEX_BNW:
    case stencilDirection::VERTEX_BW:
    case stencilDirection::VERTEX_BSW:
      return true;
    default:
      return false;
  }
}

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

inline stencilDirection makeVertexDirectionTop( const stencilDirection & dir )
{
  WALBERLA_ASSERT( isVertex( dir ) );
  switch ( dir )
  {
    case stencilDirection::VERTEX_C:
    case stencilDirection::VERTEX_BC:
      return stencilDirection::VERTEX_TC;
    case stencilDirection::VERTEX_S:
    case stencilDirection::VERTEX_BS:
      return stencilDirection::VERTEX_TS;
    case stencilDirection::VERTEX_SE:
    case stencilDirection::VERTEX_BSE:
      return stencilDirection::VERTEX_TSE;
    case stencilDirection::VERTEX_E:
    case stencilDirection::VERTEX_BE:
      return stencilDirection::VERTEX_TE;
    case stencilDirection::VERTEX_NE:
    case stencilDirection::VERTEX_BNE:
      return stencilDirection::VERTEX_TNE;
    case stencilDirection::VERTEX_N:
    case stencilDirection::VERTEX_BN:
      return stencilDirection::VERTEX_TN;
    case stencilDirection::VERTEX_NW:
    case stencilDirection::VERTEX_BNW:
      return stencilDirection::VERTEX_TNW;
    case stencilDirection::VERTEX_W:
    case stencilDirection::VERTEX_BW:
      return stencilDirection::VERTEX_TW;
    case stencilDirection::VERTEX_SW:
    case stencilDirection::VERTEX_BSW:
      return stencilDirection::VERTEX_TSW;

    default:
      return dir;
  }
}

inline stencilDirection makeVertexDirectionBottom( const stencilDirection & dir )
{
  WALBERLA_ASSERT( isVertex( dir ) );
  switch ( dir )
  {
    case stencilDirection::VERTEX_C:
    case stencilDirection::VERTEX_TC:
      return stencilDirection::VERTEX_BC;
    case stencilDirection::VERTEX_S:
    case stencilDirection::VERTEX_TS:
      return stencilDirection::VERTEX_BS;
    case stencilDirection::VERTEX_SE:
    case stencilDirection::VERTEX_TSE:
      return stencilDirection::VERTEX_BSE;
    case stencilDirection::VERTEX_E:
    case stencilDirection::VERTEX_TE:
      return stencilDirection::VERTEX_BE;
    case stencilDirection::VERTEX_NE:
    case stencilDirection::VERTEX_TNE:
      return stencilDirection::VERTEX_BNE;
    case stencilDirection::VERTEX_N:
    case stencilDirection::VERTEX_TN:
      return stencilDirection::VERTEX_BN;
    case stencilDirection::VERTEX_NW:
    case stencilDirection::VERTEX_TNW:
      return stencilDirection::VERTEX_BNW;
    case stencilDirection::VERTEX_W:
    case stencilDirection::VERTEX_TW:
      return stencilDirection::VERTEX_BW;
    case stencilDirection::VERTEX_SW:
    case stencilDirection::VERTEX_TSW:
      return stencilDirection::VERTEX_BSW;

    default:
      return dir;
  }
}

}//namespace hhg

#pragma once

#include <cassert>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"

namespace hhg {
namespace BubbleEdge {
//FIXME this can be removed after we moved into walberla namespace
using namespace walberla::mpistubs;

constexpr inline uint_t indexEdgeStencil( const stencilDirection dir )
{
   typedef hhg::stencilDirection sD;
   switch( dir )
   {
   case sD::CELL_GRAY_SW:
      return 0;
   case sD::CELL_BLUE_SE:
      return 1;
   case sD::CELL_GRAY_SE:
      return 2;
   case sD::CELL_GRAY_NW:
      return 3;
   case sD::CELL_BLUE_NW:
      return 4;
   case sD::CELL_GRAY_NE:
      return 5;
   default:
      return std::numeric_limits< size_t >::max();
   }
}

constexpr std::array< stencilDirection, 6 > neighbors = {{stencilDirection::CELL_GRAY_SE,
                                                          stencilDirection::CELL_GRAY_NE,
                                                          stencilDirection::CELL_GRAY_NW,
                                                          stencilDirection::CELL_GRAY_SW,
                                                          stencilDirection::CELL_BLUE_SE,
                                                          stencilDirection::CELL_BLUE_NW}};

constexpr std::array< stencilDirection, 3 > neighbors_south = {
    {stencilDirection::CELL_GRAY_SW, stencilDirection::CELL_BLUE_SE, stencilDirection::CELL_GRAY_SE}};

constexpr std::array< stencilDirection, 3 > neighbors_north = {
    {stencilDirection::CELL_GRAY_NW, stencilDirection::CELL_BLUE_NW, stencilDirection::CELL_GRAY_NE}};

//first face is south face by convention

constexpr inline size_t indexFaceFromVertex( const uint_t& level, size_t pos, stencilDirection dir )
{
   typedef stencilDirection sD;
   const size_t             vertexOnEdge = levelinfo::num_microvertices_per_edge( level );
   assert( pos >= 0 );
   assert( pos <= vertexOnEdge - 1 );
   const size_t startFaceS = 0;
   const size_t startFaceN = 2 * ( vertexOnEdge - 1 ) - 1;
   switch( dir )
   {
   case sD::CELL_GRAY_SE:
      return startFaceS + pos * 2;
   case sD::CELL_GRAY_NE:
      return startFaceN + pos * 2;
   case sD::CELL_GRAY_NW:
      return startFaceN + pos * 2 - 2;
   case sD::CELL_GRAY_SW:
      return startFaceS + ( pos - 1 ) * 2;
   case sD::CELL_BLUE_SE:
      return startFaceS + pos * 2 - 1;
   case sD::CELL_BLUE_NW:
      return startFaceN + pos * 2 - 1;
   default:
      // assert(false);
      return std::numeric_limits< size_t >::max();
   }
}

} // namespace BubbleEdge
} // namespace hhg

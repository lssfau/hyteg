#pragma once

#include "tinyhhg_core/macros.hpp"

namespace hhg{
namespace BubbleEdge{
//FIXME this can be removed after we moved into walberla namespace
using namespace walberla::mpistubs;

namespace EdgeCoordsVertex {
enum DirVertex {
  CELL_GRAY_SW = 0,
  CELL_BLUE_SE = 1,
  CELL_GRAY_SE = 2,
  CELL_GRAY_NW = 3,
  CELL_BLUE_NW = 4,
  CELL_GRAY_NE = 5
};

const DirVertex neighbors[] =
    {CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW, CELL_GRAY_SW,
     CELL_BLUE_SE, CELL_BLUE_NW};

const DirVertex neighbors_south[] =
    {CELL_GRAY_SW, CELL_BLUE_SE, CELL_GRAY_SE};

const DirVertex neighbors_north[] =
    {CELL_GRAY_NW, CELL_BLUE_NW, CELL_GRAY_NE};

//first face is south face by convention

template<size_t Level>
inline size_t index(size_t pos, DirVertex dir) {
  const size_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  WALBERLA_ASSERT_LESS_EQUAL(pos,vertexOnEdge);
  const size_t startFaceS = 0;
  const size_t startFaceN = 2 * (vertexOnEdge - 1) - 1;
  switch (dir) {
    case CELL_GRAY_SE:
      return startFaceS + pos * 2;
    case CELL_GRAY_NE:
      return startFaceN + pos * 2;
    case CELL_GRAY_NW:
      return startFaceN + pos * 2 - 2;
    case CELL_GRAY_SW:
      return startFaceS + (pos -1) * 2;
    case CELL_BLUE_SE:
      return startFaceS + pos * 2 -1;
    case CELL_BLUE_NW:
      return startFaceN + pos * 2 - 1;
    default:
      WALBERLA_ASSERT(false, "wrong dir");
      return std::numeric_limits<size_t>::max();
  }

}

SPECIALIZE(size_t, index, edge_index)

}//namespace EdgeCoordsVertex
}//namespace BubbleToP1Edge
}//namespace hhg

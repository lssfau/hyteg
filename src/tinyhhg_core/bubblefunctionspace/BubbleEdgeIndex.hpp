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
  WALBERLA_ABORT("Implement me.");
  return std::numeric_limits<size_t>::max();
}

SPECIALIZE(size_t, index, edge_index)

}//namespace EdgeCoordsVertex
}//namespace BubbleToP1Edge
}//namespace hhg

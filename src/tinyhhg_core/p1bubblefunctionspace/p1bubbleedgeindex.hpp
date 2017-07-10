#pragma once

#include "tinyhhg_core/macros.hpp"

namespace hhg{
namespace P1BubbleEdge{
//FIXME this can be removed after we moved into walberla namespace
using namespace walberla::mpistubs;

namespace EdgeCoordsVertex {
enum DirVertex {
  VERTEX_S  = 0,
  VERTEX_SE = 1,
  VERTEX_W  = 2,
  VERTEX_C  = 3,
  VERTEX_E  = 4,
  VERTEX_NW = 5,
  VERTEX_N  = 6,
  CELL_GRAY_SW = 7,
  CELL_BLUE_SE = 8,
  CELL_GRAY_SE = 9,
  CELL_GRAY_NW = 10,
  CELL_BLUE_NW = 11,
  CELL_GRAY_NE = 12
};

const DirVertex neighbors_with_center[13] =
    {VERTEX_C,
     VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW, CELL_GRAY_SW,
     CELL_BLUE_SE, CELL_BLUE_NW};
const DirVertex neighbors[12] =
    {VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW, CELL_GRAY_SW,
     CELL_BLUE_SE, CELL_BLUE_NW};

const DirVertex neighbors_edge[] =
    {VERTEX_E, VERTEX_W};

const DirVertex neighbors_south[] =
    {VERTEX_S, VERTEX_SE, CELL_GRAY_SW, CELL_BLUE_SE, CELL_GRAY_SE};

const DirVertex neighbors_north[] =
    {VERTEX_NW, VERTEX_N, CELL_GRAY_NW, CELL_BLUE_NW, CELL_GRAY_NE};

template<size_t Level>
inline size_t index(size_t pos, DirVertex dir) {
  const size_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  //const size_t grayBaseLength = vertexOnEdge -1;
  //const size_t blueBaseLength = vertexOnEdge -2;
  const size_t startFaceS = vertexOnEdge;
  const size_t startFaceN = 4 * (vertexOnEdge - 1);
  //const size_t totalCellGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = pos;
  //const size_t cellGrayNE = center + totalVertices - row;
  //const size_t cellBlueNW = cellGrayNE + (totalCellGray - row) -1;
  switch (dir) {
    case VERTEX_C:
      return center;
    case VERTEX_S:
      return startFaceS + pos - 1;
    case VERTEX_SE:
      return startFaceS + pos;
    case VERTEX_E:
      return center + 1;
    case VERTEX_N:
      return startFaceN + pos;
    case VERTEX_NW:
      return startFaceN + pos - 1;
    case VERTEX_W:
      return center - 1;
    case CELL_BLUE_SE:
      return startFaceS + (vertexOnEdge -1) + pos * 2 -1;
    case CELL_GRAY_NE:
      return startFaceN + (vertexOnEdge -1) + pos * 2;
    case CELL_GRAY_NW:
      return startFaceN + (vertexOnEdge -1) + pos * 2 - 2;
    case CELL_GRAY_SE:
      return startFaceS + (vertexOnEdge -1) + pos * 2;
    case CELL_BLUE_NW:
      return startFaceN + (vertexOnEdge -1) + pos * 2 - 1;
    case CELL_GRAY_SW:
      return startFaceS + (vertexOnEdge -1) + (pos -1) * 2;
  }
  return std::numeric_limits<size_t>::max();
}

SPECIALIZE(size_t, index, edge_index)

}//namespace EdgeCoordsVertex
}//namespace P1BubbleEdge
}//namespace hhg

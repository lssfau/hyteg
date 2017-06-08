#pragma once

namespace hhg
{
namespace P1BubbleFace
{
namespace CoordsVertex {
enum DirVertex {
    VERTEX_S  = 0,
    VERTEX_SE = 1,
    VERTEX_W  = 2,
    VERTEX_C  = 3,
    VERTEX_E  = 4,
    VERTEX_NW = 5,
    VERTEX_N  = 6,
    CELL_GRAY_SE = 7,
    CELL_GRAY_NW = 8,
    CELL_GRAY_NE = 9,
    CELL_BLUE_SW = 10,
    CELL_BLUE_SE = 11,
    CELL_BLUE_NW = 12
};

const DirVertex neighbors_with_center[13] =
    {VERTEX_C,
     VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW,
     CELL_BLUE_SE, CELL_BLUE_NW, CELL_BLUE_SW};
const DirVertex neighbors[12] =
    {VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW,
     CELL_BLUE_SE, CELL_BLUE_NW, CELL_BLUE_SW};

template<size_t Level>
inline size_t index(size_t row, size_t col, DirVertex dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t blueBaseLength = vertexBaseLength -2;
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t totalCellGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = (totalVertices - (vertexBaseLength-row)*(vertexBaseLength-row+1)/2) + col;
  const size_t cellGrayNE = center + totalVertices - row;
  const size_t cellBlueNW = cellGrayNE + (totalCellGray - row) -1;
  switch (dir) {
    case VERTEX_C:
      return center;
    case VERTEX_N:
      return center + vertexBaseLength - row;
    case VERTEX_E:
      return center + 1;
    case VERTEX_S:
      return center - vertexBaseLength - 1 + row;
    case VERTEX_W:
      return center - 1;
    case VERTEX_SE:
      return center - vertexBaseLength + row;
    case VERTEX_NW:
      return center + vertexBaseLength - row - 1;
    case CELL_GRAY_SE:
      return cellGrayNE - (grayBaseLength - row) -1;
    case CELL_GRAY_NE:
      return cellGrayNE;
    case CELL_GRAY_NW:
      return cellGrayNE - 1;
    case CELL_BLUE_SE:
      return cellBlueNW - (blueBaseLength - row);
    case CELL_BLUE_NW:
      return cellBlueNW;
    case CELL_BLUE_SW:
      return cellBlueNW - (blueBaseLength - row) -1;
  }
  return std::numeric_limits<size_t>::max();
}
}//namespace CoordsVertex

namespace CoordsCellGray {
enum Dir {
    VERTEX_SW = 0,
    VERTEX_SE = 1,
    VERTEX_NW = 2,
    CELL_GRAY_C = 3
};

const Dir neighbors[3] = {VERTEX_SE,VERTEX_NW,VERTEX_SW};
const Dir neighbors_with_center[4] = {CELL_GRAY_C,VERTEX_SE,VERTEX_NW,VERTEX_SW};


template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = totalVertices + totalGray - (grayBaseLength - row) * (grayBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_GRAY_C:
      return center;
    case VERTEX_SE:
      return center - totalVertices + row + 1;
    case VERTEX_SW:
      return center - totalVertices + row;
    case VERTEX_NW:
      return center - totalVertices + row + vertexBaseLength - row;
  }
  return std::numeric_limits<size_t>::max();
}


}//namesapce CoordsCellGray

namespace CoordsCellBlue {
enum Dir {
    VERTEX_SE = 0,
    VERTEX_NW = 1,
    VERTEX_NE = 2,
    CELL_BLUE_C = 3
};

const Dir neighbors[3] = {VERTEX_SE,VERTEX_NE,VERTEX_NW};
const Dir neighbors_with_center[4] = {CELL_BLUE_C,VERTEX_SE,VERTEX_NE,VERTEX_NW};


template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t blueBaseLength = vertexBaseLength -2;
  const size_t totalBlue = blueBaseLength * (blueBaseLength + 1) / 2;
  const size_t center = totalVertices + totalGray + totalBlue - (blueBaseLength - row) * (blueBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_BLUE_C:
      return center;
    case VERTEX_SE:
      return center - totalVertices + row - totalGray + row + 1;
    case VERTEX_NE:
      return center - totalVertices + row - totalGray + row + vertexBaseLength - row + 1;
    case VERTEX_NW:
      return center - totalVertices + row - totalGray + row + vertexBaseLength - row;
  }
  return std::numeric_limits<size_t>::max();
}

}//namespace CoordsCellBlue

}
}
#pragma once
#include <iterator>

namespace hhg {
namespace BubbleToP1Face {

using walberla::uint_t;

namespace CoordsVertex {
enum DirVertex {
  CELL_GRAY_SE = 0,
  CELL_GRAY_NW = 1,
  CELL_GRAY_NE = 2,
  CELL_BLUE_SW = 3,
  CELL_BLUE_SE = 4,
  CELL_BLUE_NW = 5
};

const DirVertex neighbors[] =
    {CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW,
     CELL_BLUE_SE, CELL_BLUE_NW, CELL_BLUE_SW};

template<size_t Level>
inline size_t index(size_t row, size_t col, DirVertex dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t grayBaseLength = vertexBaseLength - 1;
  const size_t blueBaseLength = vertexBaseLength - 2;
  const size_t totalVertices = vertexBaseLength*(vertexBaseLength + 1)/2;
  const size_t totalCellGray = grayBaseLength*(grayBaseLength + 1)/2;
  const size_t center = (totalVertices - (vertexBaseLength - row)*(vertexBaseLength - row + 1)/2) + col;
  const size_t cellGrayNE = center + totalVertices - row;
  const size_t cellBlueNW = cellGrayNE + (totalCellGray - row) - 1;
  switch (dir) {
    case CELL_GRAY_SE:return cellGrayNE - (grayBaseLength - row) - 1;
    case CELL_GRAY_NE:return cellGrayNE;
    case CELL_GRAY_NW:return cellGrayNE - 1;
    case CELL_BLUE_SE:return cellBlueNW - (blueBaseLength - row);
    case CELL_BLUE_NW:return cellBlueNW;
    case CELL_BLUE_SW:return cellBlueNW - (blueBaseLength - row) - 1;
  }
  return std::numeric_limits<size_t>::max();
}
}//namespace CoordsVertex

}
}
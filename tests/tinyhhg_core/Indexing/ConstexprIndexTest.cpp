
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/facedofspace/FaceDoFIndexing.hpp"

typedef size_t uint_t;
using namespace hhg;

constexpr size_t sumIndicesFace(const uint_t x, const uint_t y){
  uint_t sum = 0;
  for(uint_t i = 0; i < vertexdof::macroface::neighborsWithCenter.size(); ++i)
  {
    sum += vertexdof::macroface::indexFromVertex( 3, x, y,
                                                  vertexdof::macroface::neighborsWithCenter[i] );
  }
  return sum;
}

constexpr size_t sumBubbleFaceIndices(const uint_t x, const uint_t y){
  uint_t sum = 0;
  for(uint_t i = 0; i < facedof::macroface::neighbors.size(); ++i)
  {
    sum += facedof::macroface::indexFaceFromVertex( 3, x, y, facedof::macroface::neighbors[i] );
  }
  return sum;
}

constexpr size_t sumBubbleEdgeIndices(const uint_t x){
  uint_t sum = 0;
  for(uint_t i = 0; i < facedof::macroedge::neighbors.size(); ++i)
  {
    sum += facedof::macroedge::indexFaceFromVertex( 3, x, facedof::macroedge::neighbors[i]);
  }
  return sum;
}

constexpr size_t sumGrayDGFaceIndices(const uint_t x, const uint_t y){
  uint_t sum = 0;
  for(uint_t i = 0; i < facedof::macroface::grayFaceNeighbors.size(); ++i)
  {
    sum += facedof::macroface::indexFaceFromGrayFace( 3, x, y, facedof::macroface::grayFaceNeighbors[i] );
  }
  return sum;
}

constexpr size_t sumBlueDGFaceIndices(const uint_t x, const uint_t y){
  uint_t sum = 0;
  for(uint_t i = 0; i < facedof::macroface::blueFaceNeighbors.size(); ++i)
  {
    sum += facedof::macroface::indexFaceFromBlueFace( 3, x, y, facedof::macroface::blueFaceNeighbors[i] );
  }
  return sum;
}

int main() {
  std::vector<size_t> refOneOne = {10, 1, 2, 11, 18, 17, 9};

  static_assert( facedof::macroface::indexFaceFromVertex( 3, 1, 1, stencilDirection::CELL_GRAY_SE ) == 1,
                 "facedof::macroface::indexFaceFromVertex() failed");
  static_assert(sumBubbleFaceIndices(4,2)==194,"BubbleEdge sum failed");

  static_assert( facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_GRAY_SE ) == 8,
                 "facedof::macroedge::indexFaceFromVertex() failed");
  static_assert(sumBubbleEdgeIndices(4)==87,"BubbleEdge sum failed");

  static_assert(
  vertexdof::macroface::indexFromVertex( 3, 1, 1, stencilDirection::VERTEX_C ) == 10, "P1Face Index failed");
  static_assert(sumIndicesFace(1, 1)==68,"P1Face Index sum failed");

  static_assert( facedof::macroface::indexFaceFromGrayFace( 3, 2, 3, stencilDirection::CELL_BLUE_S ) == 51,
                 "facedof::macroface::indexFaceFromGrayFace() failed");
  static_assert(sumGrayDGFaceIndices(2, 4)==175,"P1Face Index sum failed");

  static_assert( facedof::macroface::indexFaceFromBlueFace( 3, 5, 0, stencilDirection::CELL_GRAY_N ) == 13,
                 "facedof::macroface::indexFaceFromBlueFace() failed");
  static_assert(sumBlueDGFaceIndices(6, 0)==27,"P1Face Index sum failed");

}

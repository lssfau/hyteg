
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"

typedef size_t uint_t;

constexpr size_t sumIndicesFace(const uint_t x, const uint_t y){
  uint_t sum = 0;
  for(uint_t i = 0; i < hhg::vertexdof::macroface::neighborsWithCenter.size(); ++i)
  {
    sum += hhg::vertexdof::macroface::indexFromVertex<3>(x, y, hhg::vertexdof::macroface::neighborsWithCenter[i]);
  }
  return sum;
}

constexpr size_t sumIndicesEdge(const uint_t x){
  uint_t sum = 0;
  for(uint_t i = 0; i < hhg::vertexdof::macroedge::neighborsWithCenter.size(); ++i)
  {
    sum += hhg::vertexdof::macroedge::indexFromVertex<3>(x, hhg::vertexdof::macroedge::neighborsWithCenter[i]);
  }
  return sum;
}

constexpr size_t sumBubbleFaceIndices(const uint_t x, const uint_t y){
  using namespace hhg::BubbleFace;
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors.size(); ++i)
  {
    sum += indexFaceFromVertex<3>(x, y, neighbors[i]);
  }
  return sum;
}

constexpr size_t sumBubbleEdgeIndices(const uint_t x){
  using namespace hhg::BubbleEdge;
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors.size(); ++i)
  {
    sum += indexFaceFromVertex<3>(x, neighbors[i]);
  }
  return sum;
}

constexpr size_t sumGrayDGFaceIndices(const uint_t x, const uint_t y){
  using namespace hhg::DGFace;
  uint_t sum = 0;
  for(uint_t i = 0; i < grayDGfaceneighbors.size(); ++i)
  {
    sum += indexDGFaceFromGrayDGface<3>(x, y, grayDGfaceneighbors[i]);
  }
  return sum;
}

constexpr size_t sumBlueDGFaceIndices(const uint_t x, const uint_t y){
  using namespace hhg::DGFace;
  uint_t sum = 0;
  for(uint_t i = 0; i < blueDGfaceneighbors.size(); ++i)
  {
    sum += indexDGFaceFromBlueDGface<3>(x, y, blueDGfaceneighbors[i]);
  }
  return sum;
}

int main() {
  std::vector<size_t> refOneOne = {10, 1, 2, 11, 18, 17, 9};

  static_assert(hhg::BubbleFace::indexFaceFromVertex<3>(1, 1, hhg::stencilDirection::CELL_GRAY_SE)==1,"BubbleFace Index failed");
  static_assert(sumBubbleFaceIndices(4,2)==194,"BubbleEdge sum failed");

  static_assert(hhg::BubbleEdge::indexFaceFromVertex<3>(4, hhg::stencilDirection::CELL_GRAY_SE)==8,"BubbleEdge Index failed");
  static_assert(sumBubbleEdgeIndices(4)==87,"BubbleEdge sum failed");

  static_assert(hhg::vertexdof::macroedge::indexFromVertex<3>(4, hhg::stencilDirection::VERTEX_SE)==13,"P1Edge Index failed");
  static_assert(sumIndicesEdge(3)==71,"P1Edge Index sum failed");

  static_assert(hhg::vertexdof::macroface::indexFromVertex<3>(1,1,hhg::stencilDirection::VERTEX_C)==10,"P1Face Index failed");
  static_assert(sumIndicesFace(1, 1)==68,"P1Face Index sum failed");

  static_assert(hhg::DGFace::indexDGFaceFromGrayDGface<3>(2, 3, hhg::stencilDirection::CELL_BLUE_S)==51, "DGFace Index failed");
  static_assert(sumGrayDGFaceIndices(2, 4)==175,"P1Face Index sum failed");

  static_assert(hhg::DGFace::indexDGFaceFromBlueDGface<3>(5, 0, hhg::stencilDirection::CELL_GRAY_N)==13, "DGFace Index failed");
  static_assert(sumBlueDGFaceIndices(6, 0)==27,"P1Face Index sum failed");

}

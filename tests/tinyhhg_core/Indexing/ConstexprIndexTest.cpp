#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"
#include "tinyhhg_core/dgfunctionspace/DgFaceIndex.hpp"

typedef size_t uint_t;

constexpr size_t sumIndicesFace(const uint_t x, const uint_t y){
  using namespace hhg::P1Face::FaceCoordsVertex;
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors_with_center.size(); ++i)
  {
    sum += index<3>(x, y, neighbors_with_center[i]);
  }
  return sum;
}

constexpr size_t sumIndicesEdge(const uint_t x){
  using namespace hhg::P1Edge::EdgeCoordsVertex;
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors_with_center.size(); ++i)
  {
    sum += index<3>(x, neighbors_with_center[i]);
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
  using namespace hhg::BubbleEdge::EdgeCoordsVertex;
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors.size(); ++i)
  {
    sum += index<3>(x, neighbors[i]);
  }
  return sum;
}

//constexpr uint_t sumDGFaceIndices(const uint_t x, const uint_t y){
//  return 5;
//}


int main() {
  std::vector<size_t> refOneOne = {10, 1, 2, 11, 18, 17, 9};

  static_assert(hhg::BubbleFace::indexFaceFromVertex<3>(1, 1, hhg::stencilDirection::CELL_GRAY_SE)==1,"BubbleFace Index failed");
  static_assert(sumBubbleFaceIndices(4,2)==194,"BubbleEdge sum failed");

  static_assert(hhg::BubbleEdge::EdgeCoordsVertex::index<3>(4,hhg::BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SE)==8,"BubbleEdge Index failed");
  static_assert(sumBubbleEdgeIndices(4)==87,"BubbleEdge sum failed");

  static_assert(hhg::P1Edge::EdgeCoordsVertex::index<3>(4,hhg::P1Edge::EdgeCoordsVertex::VERTEX_SE)==13,"P1Edge Index failed");
  static_assert(sumIndicesEdge(3)==71,"P1Edge Index sum failed");

  static_assert(hhg::P1Face::FaceCoordsVertex::index<3>(1,1,hhg::P1Face::FaceCoordsVertex::VERTEX_C)==10,"P1Face Index failed");
  static_assert(sumIndicesFace(1, 1)==68,"P1Face Index sum failed");

  //static_assert(hhg::DgFace::indexDGcellFromGrayDGCell<3>(2,3, hhg::stencilDirection::CELL_BLUE_S)==45, "DGFace Index failed");

}

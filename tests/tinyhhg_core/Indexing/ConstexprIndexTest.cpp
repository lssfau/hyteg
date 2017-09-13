#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"

typedef size_t uint_t;

constexpr size_t sumIndicesFace(const uint_t x, const uint_t y){
  using namespace hhg::P1Face::CoordsVertex;
  uint_t sum = 0;
  for(uint i = 0; i < neighbors_with_center.size(); ++i)
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


int main() {
#ifdef NDEBUG
  std::vector<size_t> refOneOne = {10, 1, 2, 11, 18, 17, 9};



  static_assert(hhg::P1Edge::EdgeCoordsVertex::index<3>(4,hhg::P1Edge::EdgeCoordsVertex::VERTEX_SE)==13,"P1Edge Index failed");
  static_assert(sumIndicesEdge(3)==71,"foo");

  static_assert(hhg::P1Face::CoordsVertex::index<3>(1,1,hhg::P1Face::CoordsVertex::VERTEX_C)==10,"P1Face Index failed");
  static_assert(sumIndicesFace(1, 1)==68,"foo");
#endif
}

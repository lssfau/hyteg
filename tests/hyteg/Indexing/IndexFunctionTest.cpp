#include <cstddef>
#include <array>
#include <cmath>
#include <iostream>

typedef size_t uint_t;

constexpr inline uint_t num_microvertices_per_edge(uint_t level)
{
  return (uint_t) std::pow(2, level) + 1;
}

enum DirVertex {
  VERTEX_S  = 0,
  VERTEX_SE = 1,
  VERTEX_W  = 2,
  VERTEX_C  = 3,
  VERTEX_E  = 4,
  VERTEX_NW = 5,
  VERTEX_N  = 6
};

constexpr std::array<DirVertex,7> neighbors_with_center =
  {VERTEX_C,
   VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W};

template<uint_t Level>
constexpr inline uint_t index(uint_t pos, DirVertex dir) {
  const uint_t vertexOnEdge = num_microvertices_per_edge(Level);
  const uint_t startFaceS = vertexOnEdge;
  const uint_t startFaceN = vertexOnEdge + vertexOnEdge - 1;
  const uint_t center = pos;
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
  }
}

int main(){
  static_assert(index<3>(1,VERTEX_C)==1,"foo");
  uint_t sum = 0;
  for(uint_t i = 0; i < neighbors_with_center.size(); ++i) {
    sum += index<3>(1, neighbors_with_center[0]);
  }

  for(auto neighbor : neighbors_with_center){
    sum += index<3>(1, neighbor);
  }

  printf("%zu",sum);


}
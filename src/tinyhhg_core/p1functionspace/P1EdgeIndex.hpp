#pragma once

#include <core/Abort.h>
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/levelinfo.hpp"

#include "core/DataTypes.h"
#include "core/debug/all.h"

using walberla::uint_t;

namespace hhg{
namespace P1Edge{

namespace EdgeCoordsVertex {
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
  {{VERTEX_C,
    VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};

constexpr std::array<DirVertex,6> neighbors =
  {{VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};

constexpr std::array<DirVertex,2> neighbors_on_edge =
  {{VERTEX_E, VERTEX_W}};

constexpr std::array<DirVertex,2> neighbors_south =
  {{VERTEX_S, VERTEX_SE}};

constexpr std::array<DirVertex,2> neighbors_north =
  {{VERTEX_N, VERTEX_NW}};

//first face is south face by convention

template<uint_t Level>
constexpr inline uint_t index(uint_t pos, DirVertex dir) {
  const uint_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  //the check can be reinserted if walberla supports constexpr
  //WALBERLA_ASSERT_LESS_EQUAL(pos,vertexOnEdge);
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
  //the check can be reinserted if walberla supports constexpr
  //WALBERLA_ASSERT(false, "wrong dir");
  return std::numeric_limits<uint_t>::max();
}

SPECIALIZE(uint_t, index, edge_index)

}//namespace EdgeCoordsVertex
}//namespace P1Edge
}//namespace hhg

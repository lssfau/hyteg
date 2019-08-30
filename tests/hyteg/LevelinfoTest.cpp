#include <core/debug/all.h>
#include "hyteg/Levelinfo.hpp"

using namespace hyteg;

using walberla::uint_t;

int main() {

  const uint_t edge_per_edge[] = { 1, 2, 4 , 8, 16};
  const uint_t edge_per_face[] = { 3, 9, 30, 108, 408 };
  const uint_t edge_per_cell[] = { 6, 25, 130, 804, 5576 };
  const uint_t face_per_face[] = { 1, 4, 16, 64, 256};
  const uint_t vertices_per_face[] = { 3, 6, 15, 45, 153};
  const uint_t vertices_per_edge[] = { 2, 3, 5, 9, 17};

  for (uint_t level = 0; level <= 4; ++level){
    WALBERLA_CHECK_EQUAL(levelinfo::num_microedges_per_edge(level),edge_per_edge[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microedges_per_face(level),edge_per_face[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microedges_per_cell(level),edge_per_cell[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microfaces_per_face(level),face_per_face[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microvertices_per_face(level),vertices_per_face[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microvertices_per_edge(level),vertices_per_edge[level],"level was " << level);
    WALBERLA_CHECK_EQUAL(levelinfo::num_microvertices_per_vertex(level),1,"level was " << level);
  }
}

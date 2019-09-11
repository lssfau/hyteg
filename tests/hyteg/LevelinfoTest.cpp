/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
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

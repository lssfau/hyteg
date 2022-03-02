/*
 * Copyright (c) 2022 Benjamin Mann
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

#include "simplexData.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply loadbalancing directly on our datastructures */
void loadbalancing( std::vector< uint_t >&   vertices_targetRank,
                    std::vector< EdgeData >& edges,
                    std::vector< FaceData >& faces,
                    std::vector< CellData >& cells,
                    const uint_t&            n_processes )
{
   // roundrobin
   uint_t cur = 0;

   for ( auto& rnk : vertices_targetRank )
   {
      rnk = cur % n_processes;
      ++cur;
   }
   for ( auto& edge : edges )
   {
      edge.setTargetRank( cur % n_processes );
      ++cur;
   }
   for ( auto& face : faces )
   {
      face.setTargetRank( cur % n_processes );
      ++cur;
   }
   for ( auto& cell : cells )
   {
      cell.setTargetRank( cur % n_processes );
      ++cur;
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg

/*
 * Copyright (c) 2021 Benjamin Mann
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

#pragma once

#include "simplex.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply edge bisection
      @param edge the edge to be bisected
      @param vtx index of the vertex on the edge midpoint
   */
inline void bisect_edge( std::shared_ptr< Simplex1 > edge, int64_t vtx )
{
   WALBERLA_ASSERT( not edge->has_children() );
   edge->set_midpoint_idx( vtx );
   edge->add_child( std::make_shared< Simplex1 >( edge->get_vertices()[0], vtx, edge ) );
   edge->add_child( std::make_shared< Simplex1 >( vtx, edge->get_vertices()[1], edge ) );
}

} // namespace adaptiveRefinement
} // namespace hyteg

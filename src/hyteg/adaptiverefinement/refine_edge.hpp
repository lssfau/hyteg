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
   @param coords     global coordinates of all vertices in the mesh
   @param vtxs       geometrymap, boundaryflag etc. of all vertices in the mesh
   @param edge       the edge to be bisected
   @return vtx index of the newly created vertex on the edge midpoint
*/
inline uint_t
    bisect_edge( EnumeratedList< Point3D >& coords, EnumeratedList< VertexData >& vtxs, std::shared_ptr< Simplex1 > edge )
{
   WALBERLA_ASSERT( !edge->has_children() );

   // create child IDs
   auto childIDs = edge->getPrimitiveID().createChildren();

   // vertex indices
   const auto vtx0  = edge->get_vertices()[0];
   const auto vtx1  = edge->get_vertices()[1];
   uint_t     vtx01 = vtxs.get_next_idx();

   // add midpoint to list of all vertices
   coords.append( ( coords[vtx0] + coords[vtx1] ) / 2 );
   // add vertex data
   vtxs.append(
       VertexData( edge->getGeometryMap(), edge->getBoundaryFlag(), childIDs[0], { { vtx01 } }, edge->getTargetRank() ) );

   // refine edge
   edge->set_midpoint_idx( vtx01 );
   edge->add_child( std::make_shared< Simplex1 >( vtx0, vtx01, edge ) );
   edge->add_child( std::make_shared< Simplex1 >( vtx01, vtx1, edge ) );
   // add IDs to child edges
   edge->get_children()[0]->setPrimitiveID( childIDs[1] );
   edge->get_children()[1]->setPrimitiveID( childIDs[2] );

   return vtx01;
}

} // namespace adaptiveRefinement
} // namespace hyteg

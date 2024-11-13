/*
 * Copyright (c) 2024 Marcus Mohr.
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

#include "hyteg/geometry/IdentityMap.hpp"

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void IdentityMap::setMap( SetupPrimitiveStorage& setupStorage )
{
   auto blendingMap = std::make_shared< IdentityMap >();

   for ( auto it : setupStorage.getCells() )
   {
      Cell& cell = *it.second;
      setupStorage.setGeometryMap( cell.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getFaces() )
   {
      Face& face = *it.second;
      setupStorage.setGeometryMap( face.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getEdges() )
   {
      Edge& edge = *it.second;
      setupStorage.setGeometryMap( edge.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getVertices() )
   {
      Vertex& vertex = *it.second;
      setupStorage.setGeometryMap( vertex.getID(), blendingMap );
   }
}

} // namespace hyteg

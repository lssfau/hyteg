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

#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

std::shared_ptr< Simplex1 > SimplexFactory::make_edge( uint_t a, uint_t b, Color color )
{
   auto edge = std::make_shared< Simplex1 >( _vertices[a], _vertices[b], nullptr, color, _geometryMap, _boundaryFlag );
   this->add_edge( a, b, edge );
   return edge;
}

std::shared_ptr< Simplex2 > SimplexFactory::make_face( uint_t a, uint_t b, uint_t c )
{
   std::array< uint_t, 3 >                      V{ _vertices[a], _vertices[b], _vertices[c] };
   std::array< std::shared_ptr< Simplex1 >, 3 > E{ _edges[{ a, b }], _edges[{ b, c }], _edges[{ c, a }] };

   auto face = std::make_shared< Simplex2 >( V, E, nullptr, _geometryMap, _boundaryFlag, _targetRank );
   this->add_face( a, b, c, face );
   return face;
}

std::shared_ptr< Simplex3 > SimplexFactory::make_cell( uint_t a, uint_t b, uint_t c, uint_t d )
{
   std::array< uint_t, 4 > V{ _vertices[a], _vertices[b], _vertices[c], _vertices[d] };

   std::array< std::shared_ptr< Simplex1 >, 6 > E{
       _edges[{ a, b }], _edges[{ a, c }], _edges[{ a, d }], _edges[{ b, c }], _edges[{ b, d }], _edges[{ c, d }] };

   std::array< std::shared_ptr< Simplex2 >, 4 > F{
       _faces[{ a, b, c }], _faces[{ a, b, d }], _faces[{ a, c, d }], _faces[{ b, c, d }] };

   return std::make_shared< Simplex3 >( V, E, F, _parent, _geometryMap, _boundaryFlag, _targetRank );
}

} // namespace adaptiveRefinement
} // namespace hyteg
/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include <array>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

// Version for triangles
MeshInfo MeshInfo::singleTriangle( const Point2D& v1, const Point2D& v2, const Point2D& v3 )
{
   MeshInfo meshInfo;

   meshInfo.vertices_[0] = MeshInfo::Vertex( 0, Point3D( {v1[0], v1[1], real_c( 0 )} ), 1 );
   meshInfo.vertices_[1] = MeshInfo::Vertex( 1, Point3D( {v2[0], v2[1], real_c( 0 )} ), 1 );
   meshInfo.vertices_[2] = MeshInfo::Vertex( 2, Point3D( {v3[0], v3[1], real_c( 0 )} ), 1 );

   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{0, 1}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{0, 2}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{1, 2}} ), 1 ) );

   meshInfo.addFace( Face( std::vector< IDType >( {{0, 1, 2}} ), 0 ) );

   return meshInfo;
}

// Version for tetrahedrons
MeshInfo MeshInfo::singleTetrahedron( const std::array< Point3D, 4 >& vertices )
{
   MeshInfo meshInfo;

   for ( uint_t idx = 0; idx < 4; idx++ )
   {
      meshInfo.vertices_[idx] = MeshInfo::Vertex( idx, vertices[idx], 1 );
   }

   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{0, 1}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{0, 2}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{0, 3}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{1, 2}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{1, 3}} ), 1 ) );
   meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{2, 3}} ), 1 ) );

   meshInfo.addFace( Face( std::vector< IDType >( {{0, 1, 2}} ), 1 ) );
   meshInfo.addFace( Face( std::vector< IDType >( {{0, 1, 3}} ), 1 ) );
   meshInfo.addFace( Face( std::vector< IDType >( {{0, 2, 3}} ), 1 ) );
   meshInfo.addFace( Face( std::vector< IDType >( {{1, 2, 3}} ), 1 ) );

   auto cell                           = Cell( std::vector< IDType >( {{0, 1, 2, 3}} ), 0 );
   meshInfo.cells_[cell.getVertices()] = cell;

   return meshInfo;
}

} // namespace hyteg

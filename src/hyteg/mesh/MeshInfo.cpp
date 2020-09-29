/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include <array>
#include <vector>

namespace hyteg {

void MeshInfo::setAllMeshBoundaryFlags( const uint_t& meshBoundaryFlag )
{
   for ( auto& p : vertices_ )
   {
      p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : edges_ )
   {
      p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : faces_ )
   {
      p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : cells_ )
   {
      p.second.setBoundaryFlag( meshBoundaryFlag );
   }
}

void MeshInfo::setMeshBoundaryFlagsByVertexLocation( const uint_t&                                    meshBoundaryFlag,
                                                     const std::function< bool( const Point3D& x ) >& onBoundary,
                                                     const bool&                                      allVertices )
{
   auto cond = [allVertices, onBoundary]( const std::vector< Point3D >& coordinates ) {
      if ( allVertices )
         return std::all_of( coordinates.begin(), coordinates.end(), onBoundary );
      else
         return std::any_of( coordinates.begin(), coordinates.end(), onBoundary );
   };

   for ( auto& p : vertices_ )
   {
      if ( cond( {p.second.getCoordinates()} ) )
         p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : edges_ )
   {
      std::vector< Point3D > coordinates;
      for ( auto v : p.second.getVertices() )
      {
         coordinates.push_back( vertices_[v].getCoordinates() );
      }
      if ( cond( coordinates ) )
         p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : faces_ )
   {
      std::vector< Point3D > coordinates;
      for ( auto v : p.second.getVertices() )
      {
         coordinates.push_back( vertices_[v].getCoordinates() );
      }
      if ( cond( coordinates ) )
         p.second.setBoundaryFlag( meshBoundaryFlag );
   }

   for ( auto& p : cells_ )
   {
      std::vector< Point3D > coordinates;
      for ( auto v : p.second.getVertices() )
      {
         coordinates.push_back( vertices_[v].getCoordinates() );
      }
      if ( cond( coordinates ) )
         p.second.setBoundaryFlag( meshBoundaryFlag );
   }
}

void MeshInfo::addEdge( const Edge & edge )
{
  WALBERLA_CHECK_UNEQUAL( edge.getVertices()[0], edge.getVertices()[1], "[Mesh] Mesh contains edge with zero length." );

  std::array< IDType, 2 > sortedVertexIDs;

  if ( edge.getVertices()[0] < edge.getVertices()[1] )
  {
    sortedVertexIDs = edge.getVertices();
  }
  else if ( edge.getVertices()[1] < edge.getVertices()[0] )
  {
    sortedVertexIDs[0] = edge.getVertices()[1];
    sortedVertexIDs[1] = edge.getVertices()[0];
  }

  if (edges_.count( sortedVertexIDs ) == 0)
  {
    edges_[ sortedVertexIDs ] = Edge( sortedVertexIDs, edge.getBoundaryFlag() );
  }
}


void MeshInfo::addFace( const Face & face )
{
  std::vector< IDType > sortedVertexIDs = face.getVertices();
  std::sort( sortedVertexIDs.begin(), sortedVertexIDs.end() );

  WALBERLA_CHECK_EQUAL( std::set< IDType >( sortedVertexIDs.begin(), sortedVertexIDs.end() ).size(), sortedVertexIDs.size(),
                        "[Mesh] Mesh contains face with duplicate vertices." );

  if ( faces_.count( sortedVertexIDs ) == 0 )
  {
    faces_[ sortedVertexIDs ] = Face( sortedVertexIDs, face.getBoundaryFlag() );
  }

}


void MeshInfo::addCellAndAllEdgesAndFaces( const Cell & cell )
{
  const auto cellCoordinates = cell.getVertices();

  WALBERLA_ASSERT_EQUAL( cells_.count( cellCoordinates ), 0 );
  cells_[cellCoordinates] = cell;

  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[1] }} ), 0 ) );
  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[2] }} ), 0 ) );
  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[3] }} ), 0 ) );
  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[1], cellCoordinates[2] }} ), 0 ) );
  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[1], cellCoordinates[3] }} ), 0 ) );
  addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );

  addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[1], cellCoordinates[2] }} ), 0 ) );
  addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[1], cellCoordinates[3] }} ), 0 ) );
  addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );
  addFace( Face( std::vector< IDType >( {{ cellCoordinates[1], cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );
}

} // namespace hyteg

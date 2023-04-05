/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/types/PointND.hpp"

#include <array>
#include <vector>

namespace hyteg {

using walberla::real_t;
using walberla::real_c;

MeshInfo MeshInfo::meshFaceChain( uint_t numFaces )
{
  return MeshInfo::meshFaceChain( numFaces, real_c(numFaces / 2) + real_c(numFaces % 2), real_c(1) );
}

MeshInfo MeshInfo::meshFaceChain( uint_t numFaces, real_t width, real_t height )
{
  MeshInfo meshInfo;

  if ( numFaces == 0 )
     return meshInfo;

  const real_t faceWidth  = width / ( real_c( numFaces / 2 ) + real_c( numFaces % 2 ) );
  const real_t faceHeight = height;

  uint_t lastLowerVertex;
  uint_t lastUpperVertex;

  meshInfo.vertices_[0] = MeshInfo::Vertex( 0, Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), 0 );
  meshInfo.vertices_[1] = MeshInfo::Vertex( 1, Point3D( real_c( 0 ), faceHeight, real_c( 0 ) ), 0 );
  meshInfo.vertices_[2] = MeshInfo::Vertex( 2, Point3D( faceWidth, real_c( 0 ), real_c( 0 ) ), 0 );

  meshInfo.addEdge( Edge( std::array< IDType, 2 >( { { 0, 1 } } ), 0 ) );
  meshInfo.addEdge( Edge( std::array< IDType, 2 >( { { 0, 2 } } ), 0 ) );
  meshInfo.addEdge( Edge( std::array< IDType, 2 >( { { 1, 2 } } ), 0 ) );

  meshInfo.addFace( Face( std::vector< IDType >( { { 0, 1, 2 } } ), 0 ) );

  lastLowerVertex = 2;
  lastUpperVertex = 1;

  for ( uint_t faceIdx = 1; faceIdx < numFaces; faceIdx++ )
  {
    if ( faceIdx % 2 == 0 )
    {
       const uint_t  lowerRightVertex         = lastUpperVertex + 1;
       const Point3D lowerRightVertexPosition = Point3D( real_c( ( faceIdx / 2 ) + 1 ) * faceWidth, real_c( 0 ), real_c( 0 ) );
       meshInfo.vertices_[lowerRightVertex]   = MeshInfo::Vertex(lowerRightVertex, lowerRightVertexPosition, 0);

      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastLowerVertex, lastUpperVertex }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastLowerVertex, lowerRightVertex }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastUpperVertex, lowerRightVertex }} ), 0 ) );

      meshInfo.addFace( Face( std::vector< IDType >( {{ lastLowerVertex, lastUpperVertex, lowerRightVertex }} ), 0 ) );

      lastLowerVertex = lowerRightVertex;
    }
    else
    {
      const uint_t  upperRightVertex         = lastLowerVertex + 1;
      const Point3D upperRightVertexPosition = Point3D( real_c( ( faceIdx / 2 ) + 1 ) * faceWidth, faceHeight, real_c( 0 ) );
      meshInfo.vertices_[upperRightVertex]   = MeshInfo::Vertex(upperRightVertex, upperRightVertexPosition, 0);

      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastLowerVertex, lastUpperVertex }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastLowerVertex, upperRightVertex }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ lastUpperVertex, upperRightVertex }} ), 0 ) );

      meshInfo.addFace( Face( std::vector< IDType >( {{ lastLowerVertex, lastUpperVertex, upperRightVertex }} ), 0 ) );

      lastUpperVertex = upperRightVertex;
    }
  }

  return meshInfo;
}

}
/*
 * Copyright (c) 2020 Daniel Drzisga.
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

#include "core/debug/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/Macros.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"
#include "hyteg/primitives/Face.hpp"

using walberla::real_c;
using walberla::uint_t;

namespace hyteg {
namespace vertexdof {

namespace macroface {

template < typename ValueType >
inline void projectNormal3D( uint_t                                                      level,
                             const Face&                                                 face,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdV,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdW )
{
   if ( face.getNumNeighborCells() == 2 )
   {
      WALBERLA_ABORT( "Cannot project normals if not a boundary face" );
   }

   auto dstU = face.getData( dstIdU )->getPointer( level );
   auto dstV = face.getData( dstIdV )->getPointer( level );
   auto dstW = face.getData( dstIdW )->getPointer( level );

   Point3D  normal;
   Matrix3r projection;
   Point3D  in;
   Point3D  out;

   Point3D x;
   Point3D xPhy;

   for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
   {
      x = coordinateFromIndex( level, face, it );
      face.getGeometryMap()->evalF( x, xPhy );

      normal_function( xPhy, normal );

      projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
      projection( 0, 1 ) = -normal[0] * normal[1];
      projection( 0, 2 ) = -normal[0] * normal[2];
      projection( 1, 0 ) = projection( 0, 1 );
      projection( 1, 1 ) = 1.0 - normal[1] * normal[1];
      projection( 1, 2 ) = -normal[1] * normal[2];
      projection( 2, 0 ) = projection( 0, 2 );
      projection( 2, 1 ) = projection( 1, 2 );
      projection( 2, 2 ) = 1.0 - normal[2] * normal[2];

      const uint_t idx = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];
      in[2] = dstW[idx];

      out = projection * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
      dstW[idx] = out[2];
   }
}

} // namespace macroface

namespace macroedge {

template < typename ValueType >
inline void projectNormal2D( uint_t                                                      level,
                             const Edge&                                                 edge,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdV )
{
   if ( edge.getNumNeighborFaces() == 2 )
   {
      WALBERLA_ABORT( "Cannot project normals if not a boundary edge" );
   }

   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );

   Face* faceS = storage->getFace( edge.neighborFaces()[0] );

   Point3D  normal;
   Matrix2r projection;
   Point2D  in;
   Point2D  out;

   Point3D x  = edge.getCoordinates()[0];
   real_t  h  = 1.0 / ( walberla::real_c( rowsize - 1 ) );
   Point3D dx = h * edge.getDirection();
   x += dx;
   Point3D xPhy;

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      faceS->getGeometryMap()->evalF( x, xPhy );
      normal_function( xPhy, normal );

      projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
      projection( 0, 1 ) = -normal[0] * normal[1];
      projection( 1, 0 ) = projection( 0, 1 );
      projection( 1, 1 ) = 1.0 - normal[1] * normal[1];

      in[0] = dstU[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      in[1] = dstV[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

      out = projection * in;

      dstU[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = out[0];
      dstV[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = out[1];

      x += dx;
   }
}

template < typename ValueType >
inline void projectNormal3D( uint_t                                                      level,
                             const Edge&                                                 edge,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdV,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdW )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );
   auto dstW = edge.getData( dstIdW )->getPointer( level );

   Point3D  normal;
   Matrix3r projection;
   Point3D  in;
   Point3D  out;

   Point3D x;
   Point3D xPhy;

   for ( const auto& it : vertexdof::macroedge::Iterator( level, 1 ) )
   {
      x = coordinateFromIndex( level, edge, it );
      edge.getGeometryMap()->evalF( x, xPhy );

      normal_function( xPhy, normal );

      projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
      projection( 0, 1 ) = -normal[0] * normal[1];
      projection( 0, 2 ) = -normal[0] * normal[2];
      projection( 1, 0 ) = projection( 0, 1 );
      projection( 1, 1 ) = 1.0 - normal[1] * normal[1];
      projection( 1, 2 ) = -normal[1] * normal[2];
      projection( 2, 0 ) = projection( 0, 2 );
      projection( 2, 1 ) = projection( 1, 2 );
      projection( 2, 2 ) = 1.0 - normal[2] * normal[2];

      const uint_t idx = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];
      in[2] = dstW[idx];

      out = projection * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
      dstW[idx] = out[2];
   }
}

inline void saveProjectNormalOperator2D( uint_t                                                   level,
                                         const Edge&                                              edge,
                                         const std::shared_ptr< PrimitiveStorage >&               storage,
                                         const std::function< void( const Point3D&, Point3D& ) >& normal_function,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdU,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdV,
                                         const std::shared_ptr< SparseMatrixProxy >&              mat )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );

   Face* faceS = storage->getFace( edge.neighborFaces()[0] );

   Point3D                 normal;
   std::vector< idx_t >    in( 2 );
   std::vector< idx_t >    out( 2 );

   Point3D x  = edge.getCoordinates()[0];
   real_t  h  = 1.0 / ( walberla::real_c( rowsize - 1 ) );
   Point3D dx = h * edge.getDirection();
   x += dx;
   Point3D xPhy;

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      faceS->getGeometryMap()->evalF( x, xPhy );
      normal_function( xPhy, normal );

      const auto idxU = dstU[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      const auto idxV = dstV[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

      mat->addValue( uint_c( idxU ), uint_c( idxU ), 1.0 - normal[0] * normal[0] );
      mat->addValue( uint_c( idxU ), uint_c( idxV ), -normal[0] * normal[1] );
      mat->addValue( uint_c( idxV ), uint_c( idxU ), -normal[0] * normal[1] );
      mat->addValue( uint_c( idxV ), uint_c( idxV ), 1.0 - normal[1] * normal[1] );

      x += dx;
   }
}

} // namespace macroedge

namespace macrovertex {

template < typename ValueType >
inline void projectNormal2D( uint_t                                                        level,
                             const Vertex&                                                 vertex,
                             const std::shared_ptr< PrimitiveStorage >&                    storage,
                             const std::function< void( const Point3D&, Point3D& ) >&      normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstIdV )
{
   WALBERLA_CHECK( storage->onBoundary( vertex.getID() ) );

   auto dstU = vertex.getData( dstIdU )->getPointer( level );
   auto dstV = vertex.getData( dstIdV )->getPointer( level );

   Face* faceS = storage->getFace( vertex.neighborFaces()[0] );

   Point3D xPhy;
   faceS->getGeometryMap()->evalF( vertex.getCoordinates(), xPhy );

   Point3D normal;
   normal_function( xPhy, normal );

   Matrix2r projection;
   projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
   projection( 0, 1 ) = -normal[0] * normal[1];
   projection( 1, 0 ) = projection( 0, 1 );
   projection( 1, 1 ) = 1.0 - normal[1] * normal[1];

   Point2D in;
   in[0]       = *dstU;
   in[1]       = *dstV;
   Point2D out = projection * in;

   *dstU = out[0];
   *dstV = out[1];
}

template < typename ValueType >
inline void projectNormal3D( uint_t                                                        level,
                             const Vertex&                                                 vertex,
                             const std::shared_ptr< PrimitiveStorage >&                    storage,
                             const std::function< void( const Point3D&, Point3D& ) >&      normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstIdV,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstIdW )
{
   auto dstU = vertex.getData( dstIdU )->getPointer( level );
   auto dstV = vertex.getData( dstIdV )->getPointer( level );
   auto dstW = vertex.getData( dstIdW )->getPointer( level );

   Point3D  normal;
   Matrix3r projection;
   Point3D  in;
   Point3D  out;

   Point3D x = vertex.getCoordinates();
   Point3D xPhy;
   vertex.getGeometryMap()->evalF( x, xPhy );

   normal_function( xPhy, normal );

   projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
   projection( 0, 1 ) = -normal[0] * normal[1];
   projection( 0, 2 ) = -normal[0] * normal[2];
   projection( 1, 0 ) = projection( 0, 1 );
   projection( 1, 1 ) = 1.0 - normal[1] * normal[1];
   projection( 1, 2 ) = -normal[1] * normal[2];
   projection( 2, 0 ) = projection( 0, 2 );
   projection( 2, 1 ) = projection( 1, 2 );
   projection( 2, 2 ) = 1.0 - normal[2] * normal[2];

   in[0] = dstU[0];
   in[1] = dstV[0];
   in[2] = dstW[0];

   out = projection * in;

   dstU[0] = out[0];
   dstV[0] = out[1];
   dstW[0] = out[2];
}

inline void saveProjectNormalOperator2D( uint_t                                                    level,
                                         const Vertex&                                             vertex,
                                         const std::shared_ptr< PrimitiveStorage >&                storage,
                                         const std::function< void( const Point3D&, Point3D& ) >&  normal_function,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >& dstIdU,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >& dstIdV,
                                         const std::shared_ptr< SparseMatrixProxy >&               mat )
{
   WALBERLA_CHECK( storage->onBoundary( vertex.getID() ) );

   auto dstU = vertex.getData( dstIdU )->getPointer( level );
   auto dstV = vertex.getData( dstIdV )->getPointer( level );

   Face* faceS = storage->getFace( vertex.neighborFaces()[0] );

   Point3D xPhy;
   faceS->getGeometryMap()->evalF( vertex.getCoordinates(), xPhy );

   Point3D normal;
   normal_function( xPhy, normal );

   const auto idxU = *dstU;
   const auto idxV = *dstV;

   mat->addValue( uint_c( idxU ), uint_c( idxU ), 1.0 - normal[0] * normal[0] );
   mat->addValue( uint_c( idxU ), uint_c( idxV ), -normal[0] * normal[1] );
   mat->addValue( uint_c( idxV ), uint_c( idxU ), -normal[0] * normal[1] );
   mat->addValue( uint_c( idxV ), uint_c( idxV ), 1.0 - normal[1] * normal[1] );
}

} // namespace macrovertex

} // namespace vertexdof
} // namespace hyteg

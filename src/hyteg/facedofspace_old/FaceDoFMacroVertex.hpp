/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

namespace hyteg {
namespace facedof {
namespace macrovertex {

template < typename ValueType >
inline void enumerate( Vertex&                                                       vertex,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                       const uint_t                                                  level,
                       uint_t&                                                       num )
{
   auto dst = vertex.getData( dstId )->getPointer( level );
   //for each adjacent edge there are two DoF where the first one is owned by the vertex
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      dst[i * 2] = static_cast< ValueType >( num++ );
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                                Level,
                         Vertex&                                                                      vertex,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                vertexMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcMemoryIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& f ) >& expr,
                         const std::shared_ptr< PrimitiveStorage >&                                                    storage )
{
   auto   vertexMemory = vertex.getData( vertexMemoryId )->getPointer( Level );
   uint_t rowsize      = levelinfo::num_microvertices_per_edge( Level );

   Point3D dir1;
   Point3D dir2;
   Point3D x;
   Point3D xBlend;

   std::vector< ValueType* > srcPtr;
   for ( auto src : srcMemoryIds )
   {
      srcPtr.push_back( vertex.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcMemoryIds.size() );

   for ( auto faceIt : vertex.neighborFaces() )
   {
      Face*  face           = storage->getFace( faceIt.getID() );
      uint_t vertexIdOnFace = face->vertex_index( vertex.getID() );
      dir1 =
          ( face->getCoordinates()[( vertexIdOnFace + 1 ) % 3] - vertex.getCoordinates() ) / ( walberla::real_c( rowsize - 1 ) );
      dir2 =
          ( face->getCoordinates()[( vertexIdOnFace + 2 ) % 3] - vertex.getCoordinates() ) / ( walberla::real_c( rowsize - 1 ) );
      x = vertex.getCoordinates() + 1.0 / 3.0 * ( dir1 + dir2 );
      for ( size_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][vertex.face_index( face->getID() ) * 2];
      }

      face->getGeometryMap()->evalF( x, xBlend );
      vertexMemory[vertex.face_index( face->getID() ) * 2] = expr( xBlend, srcVector );
   }
}

template < typename ValueType >
inline void assign( const uint_t&                                                                Level,
                    Vertex&                                                                      vertex,
                    const std::vector< ValueType >&                                              scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId )
{
   auto dst = vertex.getData( dstId )->getPointer( Level );

   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      uint_t index = i * 2;
      //tmp is necessary since dstId can also be in srcIds
      ValueType tmp = scalars[0] * vertex.getData( srcIds[0] )->getPointer( Level )[index];
      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] * vertex.getData( srcIds[k] )->getPointer( Level )[index];
      }
      dst[index] = tmp;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                                Level,
                 Vertex&                                                                      vertex,
                 const std::vector< ValueType >&                                              scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId )
{
   auto dst = vertex.getData( dstId )->getPointer( Level );

   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      uint_t index = i * 2;
      //tmp is necessary since dstId can also be in srcIds
      ValueType tmp = scalars[0] * vertex.getData( srcIds[0] )->getPointer( Level )[index];
      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] * vertex.getData( srcIds[k] )->getPointer( Level )[index];
      }
      dst[index] += tmp;
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                                level,
                             Vertex&                                                                      vertex,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId )
{
   ValueType*                dstPtr = vertex.getData( dstId )->getPointer( level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( vertex.getData( src )->getPointer( level ) );
   }

   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      uint_t    index = i * 2;
      ValueType tmp   = srcPtr[0][index];
      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp *= srcPtr[k][index];
      }
      dstPtr[index] = tmp;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                 level,
                 Vertex&                                                       vertex,
                 const ValueType                                               scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId )
{
   ValueType* dstPtr = vertex.getData( dstId )->getPointer( level );
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
      dstPtr[i * 2] += scalar;
}

template < typename ValueType >
inline ValueType
    getMaxValue( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   ValueType localMax = -std::numeric_limits< ValueType >::max();

   ValueType* srcPtr = vertex.getData( srcId )->getPointer( level );
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
      localMax = std::max( localMax, srcPtr[i * 2] );

   return localMax;
}

template < typename ValueType >
inline ValueType
    getMinValue( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   ValueType localMin = +std::numeric_limits< ValueType >::max();

   ValueType* srcPtr = vertex.getData( srcId )->getPointer( level );
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
      localMin = std::min( localMin, srcPtr[i * 2] );

   return localMin;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   ValueType localMax = ValueType( 0.0 );

   ValueType* srcPtr = vertex.getData( srcId )->getPointer( level );
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
      localMax = std::max( localMax, std::abs( srcPtr[i * 2] ) );

   return localMax;
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                                 level,
                      Vertex&                                                       vertex,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& lhsMemoryId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& rhsMemoryId )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   ValueType* lhsPtr = vertex.getData( lhsMemoryId )->getPointer( level );
   ValueType* rhsPtr = vertex.getData( rhsMemoryId )->getPointer( level );
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
      scalarProduct += rhsPtr[i * 2] * lhsPtr[i * 2];

   return scalarProduct.get();
}

template < typename ValueType >
inline void swap( const uint_t&                                                 level,
                  Vertex&                                                       vertex,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstID )
{
   vertex.getData( srcID )->swap( *vertex.getData( dstID ), level );
}

} // namespace macrovertex
} // namespace facedof
} //namespace hyteg

/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/Levelinfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

#ifdef DEBUG_ELEMENTWISE
#include "hyteg/format.hpp"
#endif

namespace hyteg {
namespace vertexdof {
namespace macrovertex {

template < typename ValueType >
inline void interpolate( const uint_t&                                                 level,
                         const Vertex&                                                 vertex,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& vertexMemoryId,
                         const ValueType&                                              scalar )
{
   auto vertexMemory = vertex.getData( vertexMemoryId )->getPointer( level );
   vertexMemory[0]   = scalar;
}

template < typename ValueType >
inline void interpolate( Vertex&                                                                      vertex,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                vertexMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr,
                         uint_t                                                                                      level )
{
   FunctionMemory< ValueType >* vertexMemory = vertex.getData( vertexMemoryId );
   std::vector< ValueType >     srcVector( srcIds.size() );

   for ( uint_t k = 0; k < srcIds.size(); ++k )
   {
      srcVector[k] = vertex.getData( srcIds[k] )->getPointer( level )[0];
   }

   Point3D xBlend;
   vertex.getGeometryMap()->evalF( vertex.getCoordinates(), xBlend );
   vertexMemory->getPointer( level )[0] = expr( xBlend, srcVector );
}

template < typename ValueType >
inline void swap( const uint_t&                                                 level,
                  Vertex&                                                       vertex,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstID )
{
   auto srcData = vertex.getData( srcID );
   auto dstData = vertex.getData( dstID );
   srcData->swap( *dstData, level );
}

template < typename ValueType >
inline void assign( Vertex&                                                                      vertex,
                    const std::vector< ValueType >&                                              scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId,
                    size_t                                                                       level )
{
   ValueType tmp = scalars[0] * vertex.getData( srcIds[0] )->getPointer( level )[0];

   for ( size_t i = 1; i < srcIds.size(); ++i )
   {
      tmp += scalars[i] * vertex.getData( srcIds[i] )->getPointer( level )[0];
   }

   vertex.getData( dstId )->getPointer( level )[0] = tmp;
}

template < typename ValueType >
inline void add( const Vertex&                                                 vertex,
                 const ValueType&                                              scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                 const uint_t&                                                 level )
{
   vertex.getData( dstId )->getPointer( level )[0] += scalar;
}

template < typename ValueType >
inline void add( Vertex&                                                                      vertex,
                 const std::vector< ValueType >&                                              scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId,
                 size_t                                                                       level )
{
   auto tmp = ValueType( 0 );

   for ( size_t i = 0; i < srcIds.size(); ++i )
   {
      tmp += scalars[i] * vertex.getData( srcIds[i] )->getPointer( level )[0];
   }

   vertex.getData( dstId )->getPointer( level )[0] += tmp;
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                                level,
                             Vertex&                                                                      vertex,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >&                dstId )
{
   ValueType tmp = vertex.getData( srcIds[0] )->getPointer( level )[0];

   for ( size_t i = 1; i < srcIds.size(); ++i )
   {
      tmp *= vertex.getData( srcIds[i] )->getPointer( level )[0];
   }

   vertex.getData( dstId )->getPointer( level )[0] = tmp;
}

template < typename ValueType >
inline ValueType dot( Vertex&                                                       vertex,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& lhsMemoryId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& rhsMemoryId,
                      size_t                                                        level )
{
   return vertex.getData( lhsMemoryId )->getPointer( level )[0] * vertex.getData( rhsMemoryId )->getPointer( level )[0];
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                                 level,
                      const Vertex&                                                 vertex,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dataID,
                      const bool&                                                   absolute )
{
   if ( absolute )
   {
      return std::abs( vertex.getData( dataID )->getPointer( level )[0] );
   }
   return vertex.getData( dataID )->getPointer( level )[0];
}

template < typename ValueType >
inline void apply( Vertex&                                                       vertex,
                   const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                   size_t                                                        level,
                   UpdateType                                                    update )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto src      = vertex.getData( srcId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );

   if ( update == Replace )
   {
      dst[0] = opr_data[0] * src[0];
   }
   else if ( update == Add )
   {
      dst[0] += opr_data[0] * src[0];
   }

   for ( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] += opr_data[i + 1] * src[i + 1];
   }
}

template < typename ValueType >
inline void smooth_gs( Vertex&                                                       vertex,
                       const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& rhsId,
                       size_t                                                        level )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );
   auto rhs      = vertex.getData( rhsId )->getPointer( level );

   dst[0] = rhs[0];

   for ( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] -= opr_data[i + 1] * dst[i + 1];
   }

   dst[0] /= opr_data[0];
}

template < typename ValueType >
inline void smooth_sor( Vertex&                                                       vertex,
                        const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& rhsId,
                        size_t                                                        level,
                        ValueType                                                     relax )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );
   auto rhs      = vertex.getData( rhsId )->getPointer( level );

   ValueType tmp;
   tmp = rhs[0];

   for ( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      tmp -= opr_data[i + 1] * dst[i + 1];
   }

   dst[0] = ( 1.0 - relax ) * dst[0] + relax * tmp / opr_data[0];
}

template < typename ValueType >
inline void smooth_jac( Vertex&                                                       vertex,
                        const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& rhsId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& tmpId,
                        size_t                                                        level )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );
   auto rhs      = vertex.getData( rhsId )->getPointer( level );
   auto tmp      = vertex.getData( tmpId )->getPointer( level );

   dst[0] = rhs[0];

   for ( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] -= opr_data[i + 1] * tmp[i + 1];
   }

   dst[0] /= opr_data[0];
}

template < typename ValueType >
inline void
    enumerate( size_t level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId, ValueType& num )
{
   auto dst = vertex.getData( dstId )->getPointer( level );
   dst[0]   = num++;
}

template < typename ValueType >
inline ValueType
    getMaxValue( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   auto src = vertex.getData( srcId )->getPointer( level );
   return src[0];
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   auto src = vertex.getData( srcId )->getPointer( level );
   return std::abs( src[0] );
}

template < typename ValueType >
inline ValueType
    getMinValue( const uint_t& level, Vertex& vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId )
{
   auto src = vertex.getData( srcId )->getPointer( level );
   return src[0];
}

template < typename ValueType >
inline void saveOperator( Vertex&                                                      vertex,
                          const PrimitiveDataID< StencilMemory< ValueType >, Vertex >& operatorId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >&    srcId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >&    dstId,
                          const std::shared_ptr< SparseMatrixProxy >&                  mat,
                          uint_t                                                       level )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto src      = vertex.getData( srcId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );

   for ( uint_t i = 0; i < vertex.getNumNeighborEdges() + 1; i++ )
   {
      mat->addValue( uint_c( dst[0] ), uint_c( src[i] ), opr_data[i] );
   }
}

inline void saveIdentityOperator( Vertex&                                                   vertex,
                                  const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >& dstId,
                                  const std::shared_ptr< SparseMatrixProxy >&               mat,
                                  uint_t                                                    level )
{
   auto dst = vertex.getData( dstId )->getPointer( level );
   mat->addValue( uint_c( dst[0] ), uint_c( dst[0] ), 1.0 );
}

template < typename ValueType >
inline void createVectorFromFunction( const Vertex&                                                 vertex,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                         vec,
                                      uint_t                                                        level )
{
   auto  src       = vertex.getData( srcId )->getPointer( level );
   idx_t numerator = vertex.getData( numeratorId )->getPointer( level )[0];

   vec->setValue( uint_c( numerator ), src[0] );
}

template < typename ValueType >
inline void createFunctionFromVector( Vertex&                                                       vertex,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                         vec,
                                      uint_t                                                        level )
{
   idx_t numerator                                 = vertex.getData( numeratorId )->getPointer( level )[0];
   vertex.getData( srcId )->getPointer( level )[0] = vec->getValue( uint_c( numerator ) );
}

inline void applyDirichletBC( Vertex&                                                   vertex,
                              std::vector< idx_t >&                                     mat,
                              uint_t                                                    level,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Vertex >& numeratorId )
{
   mat.push_back( vertex.getData( numeratorId )->getPointer( level )[0] );
}

} // namespace macrovertex
} // namespace vertexdof
} // namespace hyteg

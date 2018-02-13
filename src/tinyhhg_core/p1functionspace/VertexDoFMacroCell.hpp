
#pragma once

#include "core/debug/all.h"
#include "core/DataTypes.h"

#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/types/flags.hpp"

namespace hhg {
namespace vertexdof {
namespace macrocell {

using walberla::uint_t;
using walberla::real_t;
using walberla::real_c;

using indexing::Index;

template< uint_t Level >
inline Point3D coordinateFromIndex( const Cell & cell, const Index & index )
{
  const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( Level );
  const Point3D xStep         = ( cell.getCoordinates()[1] - cell.getCoordinates()[0] ) * stepFrequency;
  const Point3D yStep         = ( cell.getCoordinates()[2] - cell.getCoordinates()[0] ) * stepFrequency;
  const Point3D zStep         = ( cell.getCoordinates()[3] - cell.getCoordinates()[0] ) * stepFrequency;
  return cell.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() ) + zStep * real_c( index.z() );
}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl( const Cell & cell,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& cellMemoryId,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                             const std::function< ValueType( const hhg::Point3D &, const std::vector< ValueType > & )> & expr)
{
  ValueType * cellData = cell.getData( cellMemoryId )->getPointer( Level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData(src)->getPointer( Level ) );
  }

  std::vector<ValueType> srcVector( srcIds.size() );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const Point3D coordinate = coordinateFromIndex< Level >( cell, it );
    const uint_t  idx        = vertexdof::macrocell::indexFromVertex<Level>( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVector[ k ] = srcPtr[ k ][ idx ];
    }
    cellData[ idx ] = expr( coordinate, srcVector );
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate);

template< typename ValueType, uint_t Level >
inline void assignTmpl( const Cell & cell,
                        const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( Level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData( src )->getPointer( Level ));
  }

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex< Level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    ValueType tmp = scalars[0] * srcPtr[0][ idx ];

    for ( uint_t k = 1; k < srcIds.size(); ++k )
    {
      tmp += scalars[k] * srcPtr[k][ idx ];
    }
    dst[ idx ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign);


template< typename ValueType, uint_t Level >
inline void addTmpl( const Cell & cell,
                     const std::vector< ValueType > & scalars,
                     const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( Level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData( src )->getPointer( Level ));
  }

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex< Level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    ValueType tmp = scalars[0] * srcPtr[0][ idx ];

    for ( uint_t k = 1; k < srcIds.size(); ++k )
    {
      tmp += scalars[k] * srcPtr[k][ idx ];
    }
    dst[ idx ] += tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, addTmpl, add);


template< typename ValueType, uint_t Level >
inline real_t dotTmpl( const Cell & cell,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & lhsId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & rhsId)
{
  real_t sp = 0.0;

  const ValueType * lhsPtr = cell.getData( lhsId )->getPointer( Level );
  const ValueType * rhsPtr = cell.getData( rhsId )->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex< Level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
    sp += lhsPtr[ idx ] * rhsPtr[ idx ];
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE(real_t, dotTmpl, dot);


template< typename ValueType, uint_t Level >
inline void applyTmpl( Cell & cell,
                        const PrimitiveDataID< StencilMemory< ValueType >,  Cell > & operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & srcId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId,
                        const UpdateType update )
{
  typedef stencilDirection sd;

  const ValueType * operatorData = cell.getData( operatorId )->getPointer( Level );
  const ValueType * src          = cell.getData( srcId )->getPointer( Level );
        ValueType * dst          = cell.getData( dstId )->getPointer( Level );

  ValueType tmp;

  if( update == Replace )
  {
    for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
    {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex< Level >( x, y, z, sd::VERTEX_C );

      tmp = operatorData[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] * src[ centerIdx ];

      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
        const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );
        const uint_t idx        = vertexdof::macrocell::indexFromVertex< Level >( x, y, z, neighbor );
        tmp += operatorData[ stencilIdx ] * src[ idx ];
      }

      dst[ centerIdx ] = tmp;
    }
  }
  else
  {
    for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
    {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex< Level >( x, y, z, sd::VERTEX_C );

      tmp = operatorData[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] * src[ centerIdx ];

      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
        const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );
        const uint_t idx        = vertexdof::macrocell::indexFromVertex< Level >( x, y, z, neighbor );
        tmp += operatorData[ stencilIdx ] * src[ idx ];
      }

      dst[ centerIdx ] += tmp;
    }
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyTmpl, apply)


} // namespace macrocell
} // namespace vertexdof
} // namespace hhg

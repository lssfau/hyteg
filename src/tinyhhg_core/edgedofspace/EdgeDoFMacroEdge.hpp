
#pragma once

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMemory.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {
namespace edgedof {
namespace macroedge {

using walberla::uint_t;
using walberla::real_c;

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Edge & edge,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                            std::function< ValueType( const hhg::Point3D &, const std::vector<ValueType>& ) > & expr)
{
  auto edgeData = edge.getData( edgeMemoryId )->getPointer( Level );

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(edge.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcIds.size());

  const Point3D leftCoords  = edge.getCoordinates()[0];
  const Point3D rightCoords = edge.getCoordinates()[1];

  const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( Level ) );

  for ( const auto & it : indexing::edgedof::macroedge::Iterator( Level ) )
  {
    const Point3D currentCoordinates = leftCoords + microEdgeOffset + 2 * it.col() * microEdgeOffset;

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVector[k] = srcPtr[k][indexing::edgedof::macroedge::horizontalIndex< Level >( it.col() )];
    }

    edgeData[ indexing::edgedof::macroedge::indexFromHorizontalEdge< Level >( it.col(), stencilDirection::EDGE_HO_C ) ] = expr( currentCoordinates, srcVector );
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate );


template< typename ValueType, uint_t Level >
inline void addTmpl( Edge & edge, const std::vector< ValueType > & scalars,
                     const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > & srcIds,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = edge.getData( dstId )->getPointer( Level );

  for ( uint_t i = 0; i < scalars.size(); i++ )
  {
    const real_t scalar  = scalars[i];
    auto         srcData = edge.getData( srcIds[i] )->getPointer( Level );

    for ( const auto & it : indexing::edgedof::macroedge::Iterator( Level ) )
    {
      const uint_t idx = indexing::edgedof::macroedge::indexFromHorizontalEdge< Level >( it.col(), stencilDirection::EDGE_HO_C );
      dstData[ idx ] += scalar * srcData[ idx ];
    }
  }
}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add );


template< typename ValueType, uint_t Level >
inline void assignTmpl( Edge & edge, const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = edge.getData( dstId )->getPointer( Level );

  for ( const auto & it : indexing::edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = indexing::edgedof::macroedge::indexFromHorizontalEdge< Level >( it.col(), stencilDirection::EDGE_HO_C );
    dstData[ idx ] = static_cast< ValueType >( 0 );
  }

  addTmpl< ValueType, Level >( edge, scalars, srcIds, dstId );
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign );


template< typename ValueType, uint_t Level >
inline real_t dotTmpl( Edge & edge,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& lhsId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId )
{
  auto lhsData = edge.getData( lhsId )->getPointer( Level );
  auto rhsData = edge.getData( rhsId )->getPointer( Level );

  real_t scalarProduct = real_c( 0 );

  for ( const auto & it : indexing::edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = indexing::edgedof::macroedge::indexFromHorizontalEdge< Level >( it.col(), stencilDirection::EDGE_HO_C );
    scalarProduct += lhsData[ idx ] * rhsData[ idx ];
  }

  return scalarProduct;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot );


template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &dstId,
                          uint_t& num)
{
  ValueType *dst = edge.getData(dstId)->getPointer(Level);

  for(uint_t i = 0 ; i < levelinfo::num_microedges_per_edge( Level ) ; ++i){
    dst[hhg::indexing::edgedof::macroedge::horizontalIndex< Level >(i)] = num;
    ++num;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate );

} ///namespace macroedge
} ///namespace edgedof
} ///namespace hhg

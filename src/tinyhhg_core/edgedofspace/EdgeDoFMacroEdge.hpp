
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
                            std::function< ValueType( const hhg::Point3D & ) > & expr)
{
  auto edgeData = edge.getData( edgeMemoryId )->getPointer( Level );

  const Point3D leftCoords  = edge.getCoordinates()[0];
  const Point3D rightCoords = edge.getCoordinates()[1];

  const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( Level ) );

  for ( const auto & idx : indexing::edgedof::macroedge::Iterator( Level ) )
  {
    const Point3D currentCoordinates = leftCoords + microEdgeOffset + 2 * idx.col() * microEdgeOffset;
    edgeData[ indexing::edgedof::macroedge::indexFromHorizontalEdge< Level >( idx.col(), stencilDirection::EDGE_HO_C ) ] = expr( currentCoordinates );
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate );

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

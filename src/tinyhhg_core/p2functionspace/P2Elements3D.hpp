
#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace P2Elements {
namespace P2Elements3D {


inline std::set< indexing::IndexIncrement > getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( const edgedof::EdgeDoFOrientation & edgeDoFOrientation )
{
  std::set< indexing::IndexIncrement > neighborEdgeDoFs;
  for ( const auto & element : P1Elements::P1Elements3D::allCellsAtInnerVertex )
  {
    WALBERLA_ASSERT_EQUAL( element[0], stencilDirection::VERTEX_C );
    WALBERLA_ASSERT_EQUAL( element.size(), 4 );
    for ( const auto & dir0 : element )
    {
      for ( const auto & dir1 : element )
      {
        if ( dir0 != dir1 )
        {
          const auto elementEdgeOrientation = edgedof::calcEdgeDoFOrientation( vertexdof::logicalIndexOffsetFromVertex( dir0 ), vertexdof::logicalIndexOffsetFromVertex( dir1 ));

          if ( edgeDoFOrientation != elementEdgeOrientation )
          {
            continue;
          }

          neighborEdgeDoFs.insert( edgedof::calcEdgeDoFIndex( vertexdof::logicalIndexOffsetFromVertex( dir0 ), vertexdof::logicalIndexOffsetFromVertex( dir1 )));
        }
      }
    }
  }
  return neighborEdgeDoFs;
}


inline std::set< indexing::IndexIncrement > getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                          const edgedof::EdgeDoFOrientation & leafOrientation )
{
  std::set< indexing::IndexIncrement > neighborEdgeDoFs;
  const auto neighborVertices = edgedof::calcNeighboringVertexDoFIndices( centerOrientation );
  const auto secondDirectionInElements = vertexdof::stencilDirectionFromLogicalOffset( neighborVertices.at(1) - neighborVertices.at(0) );

  for ( const auto & element : P1Elements::P1Elements3D::allCellsAtInnerVertex )
  {
    WALBERLA_ASSERT_EQUAL( element[0], stencilDirection::VERTEX_C );
    WALBERLA_ASSERT_EQUAL( element.size(), 4 );

    if ( std::find( element.begin(), element.end(), secondDirectionInElements ) == element.end() )
    {
      continue;
    }

    for ( const auto & dir0 : element )
    {
      for ( const auto & dir1 : element )
      {
        if ( dir0 != dir1 )
        {
          const auto elementEdgeOrientation = edgedof::calcEdgeDoFOrientation( vertexdof::logicalIndexOffsetFromVertex( dir0 ), vertexdof::logicalIndexOffsetFromVertex( dir1 ));

          if ( leafOrientation != elementEdgeOrientation )
          {
            continue;
          }

          neighborEdgeDoFs.insert( edgedof::calcEdgeDoFIndex( vertexdof::logicalIndexOffsetFromVertex( dir0 ), vertexdof::logicalIndexOffsetFromVertex( dir1 )) + neighborVertices.at(0));
        }
      }
    }
  }
  return neighborEdgeDoFs;
}


}
}
}
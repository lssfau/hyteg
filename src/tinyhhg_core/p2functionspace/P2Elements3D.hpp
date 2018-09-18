
#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace P2Elements {
namespace P2Elements3D {

/// \brief Returns the index offsets to the neighboring edge unknowns of a micro-vertex.
///
/// Given a logical index I of a micro-vertex, this function returns the
/// index offsets that have to be added to I to reach the absolute logical indices of the
/// neighboring edgedofs of the specified orientation.
///
/// \param edgeDoFOrientation the orientation of the edges that shall be reached from the vertexdof
/// \return a set that contains all index offsets to neighboring edgedofs
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


/// \brief Returns the index offsets to the neighboring edge unknowns of a micro-edge.
///
/// Given a logical index I and an orientation of the corresponding edgedof (== centerOrientation),
/// this function returns the index offsets that have to be added to I to reach the absolute logical indices of the
/// neighboring edgedofs of a certain orientation (== leafOrientation).
///
/// \param centerOrientation the orientation of the edge from which the neighbor edges shall be reached
/// \param leafOrientation the orientation of the edges that shall be reached
/// \return a set that contains all index offsets to neighboring edgedofs with leafOrientation
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


/// \brief Returns the index offsets to the neighboring vertex unknowns of a micro-edge.
///
/// Given a logical index I of a micro-edge and its orientation, this function returns the
/// index offsets that have to be added to I to reach the absolute logical indices of the
/// neighboring vertexdofs.
///
/// \param centerOrientation the orientation of the edges from which the neighbor vertices shall be reached
/// \return a set that contains all index offsets to neighboring vertexdofs
inline std::set< indexing::IndexIncrement > getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( const edgedof::EdgeDoFOrientation & centerOrientation )
{
  std::set< indexing::IndexIncrement > neighborVertexDoFs;
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

    for ( const auto & dir : element )
    {
      neighborVertexDoFs.insert( vertexdof::logicalIndexOffsetFromVertex( dir ) + neighborVertices.at(0) );
    }
  }
  return neighborVertexDoFs;
}


}
}
}
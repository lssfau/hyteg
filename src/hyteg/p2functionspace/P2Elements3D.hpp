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

#pragma once

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/forms/P2Form.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"

#include <set>

namespace hyteg {
namespace P2Elements {
namespace P2Elements3D {

/// \brief Returns the (sorted) array indices of the two micro-vertices of the given element that span the edge with the specified orientation there are any.
inline std::optional< std::array< uint_t, 2 > > edgeWithOrientationFromElement( const std::array< indexing::Index, 4 > & elementVertices,
                                                                                     const edgedof::EdgeDoFOrientation & orientation )
{
  for ( uint_t vertex0 = 0; vertex0 < 4; vertex0++ )
  {
    for ( uint_t vertex1 = 0; vertex1 < vertex0; vertex1++ )
    {
      if ( edgedof::calcEdgeDoFOrientation( elementVertices.at( vertex0 ), elementVertices.at( vertex1 ) ) == orientation )
        return std::array< uint_t, 2 >({ vertex1, vertex0 });
    }
  }
  return {};
}


/// \brief Returns the index offsets to the neighboring edge unknowns of a micro-vertex.
///
/// Given a logical index I of a micro-vertex, this function returns the
/// index offsets that have to be added to I to reach the absolute logical indices of the
/// neighboring edgedofs of the specified orientation.
///
/// \param edgeDoFOrientation the orientation of the edges that shall be reached from the vertexdof
/// \return a set that contains all index offsets to neighboring edgedofs
inline std::set< indexing::Index > getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( const edgedof::EdgeDoFOrientation & edgeDoFOrientation )
{
  std::set< indexing::Index > neighborEdgeDoFs;
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
inline std::set< indexing::Index > getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                          const edgedof::EdgeDoFOrientation & leafOrientation )
{
  std::set< indexing::Index > neighborEdgeDoFs;
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
inline std::set< indexing::Index > getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( const edgedof::EdgeDoFOrientation & centerOrientation )
{
  std::set< indexing::Index > neighborVertexDoFs;
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


/// \brief Calculates the edge to vertex stencil weights for one (leaf-)orientation from the stiffness matrices of neighboring
///        elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-vertex' location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
///
/// \param microVertexIndex the logical index of the micro-vertex in a macro-cell (can also be located on the macro-cell's boundary)
/// \param leafOrientation the orientation of the stencil leaves
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param ufcGen the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from index offsets to stencil weights
///
template< typename UFCOperator >
inline std::map< indexing::Index, real_t > calculateEdgeToVertexStencilInMacroCell( const indexing::Index & microVertexIndex,
                                                                                             const edgedof::EdgeDoFOrientation & leafOrientation,
                                                                                             const Cell & cell, const uint_t & level,
                                                                                             const UFCOperator & ufcGen )
{
  typedef stencilDirection sd;
  std::map< indexing::Index, real_t > macroCellStencilEntries;

  const auto neighboringElements = P1Elements::P1Elements3D::getNeighboringElements( microVertexIndex, level );

  // 1. Going over all neighboring cells of a micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
  //    VERTEX_C represents the reference micro-vertex.
  for ( const auto & cellAtVertex : neighboringElements )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    const std::array< indexing::Index, 4 > elementAsIndices = {
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[0] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[1] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[2] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[3] ),
    };

    // Since there are no parallel edges in the same element, there is at most one edge with the specified orientation in the current element.
    // If there is none, we don't do anything with this element.
    if ( const auto optEdge = edgeWithOrientationFromElement( elementAsIndices, leafOrientation ))
    {
      const auto edge = optEdge.value();

      // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
      std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ ) {
        logicalOffsetsFromCenter[localID] = microVertexIndex + elementAsIndices.at( localID );
      }

      // 5. Adding contribution to stencil

      // obtain micro cell coordinates
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ ) {
        geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      // find index into stencil
      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( elementAsIndices.at( edge.at(0) ), elementAsIndices.at( edge.at(1) ) );

      // obtain entry of local element matrix and add to stencil weight
      macroCellStencilEntries[ edgeDoFIndex ] += ufcGen.integrate( geometricCoordinates, { 0, 0 }, edge );

    }
  }

  return macroCellStencilEntries;
}


/// \brief Calculates the vertex to edge stencil weights for one orientation (of the center edge) from the stiffness matrices of neighboring elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-edge's location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
///
/// \param microEdgeIndex the logical index of the micro-edge in a macro-cell (can also be located on the macro-cell's boundary)
/// \param centerOrientation the orientation of the stencil center micro-edge
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param ufcGen the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from index offsets to stencil weights
///
template< typename UFCOperator >
inline std::map< indexing::Index, real_t > calculateVertexToEdgeStencilInMacroCell( const indexing::Index & microEdgeIndex,
                                                                                             const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                             const Cell & cell, const uint_t & level,
                                                                                             const UFCOperator & ufcGen )
{
  typedef stencilDirection                     sd;
  std::map< indexing::Index, real_t > macroCellStencilEntries;

  const auto  offsetsToNeighborVertices = edgedof::calcNeighboringVertexDoFIndices( centerOrientation );
  const auto& neighboringVertex0        = offsetsToNeighborVertices.at( 0 );
  const auto  neighboringElementsAtVertex0 =
      P1Elements::P1Elements3D::getNeighboringElements( microEdgeIndex + neighboringVertex0, level );
  const auto secondDirectionInElements =
      vertexdof::stencilDirectionFromLogicalOffset( offsetsToNeighborVertices.at( 1 ) - offsetsToNeighborVertices.at( 0 ) );

  // 1. Going over all neighboring cells of the first neighboring micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions.
  for ( const auto& cellAtVertex : neighboringElementsAtVertex0 )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    // Exclude elements that do not contain the second vertex that defines the current edge.
    if ( std::find( cellAtVertex.begin(), cellAtVertex.end(), secondDirectionInElements ) == cellAtVertex.end() )
    {
      continue;
    }

    const std::array< indexing::Index, 4 > elementAsIndices = {
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[0] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[1] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[2] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[3] ),
    };

    // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex.
    //    The reference micro-index is the 'first' of the two indices that define the current edge.
    std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      logicalOffsetsFromCenter[localID] = microEdgeIndex + neighboringVertex0 + elementAsIndices.at( localID );
    }

    // 5. Adding contribution to stencil

    // obtain micro cell coordinates
    std::array< Point3D, 4 > geometricCoordinates;
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
    }

    // find vertices of centerEdge
    auto centerEdge = edgeWithOrientationFromElement( elementAsIndices, centerOrientation ).value();

    if ( std::is_same< UFCOperator, P2Form >::value )
    {
      // obtain weights from local element matrix
      std::vector< P2Form::dofPosByVertexPair3D > leafPos = { { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 3 } } };
      std::vector<real_t> weights = ufcGen.integrate( geometricCoordinates, centerEdge, leafPos );

      // add values at correct index position into stencil
      for ( uint_t localVertexID = 0; localVertexID < weights.size(); localVertexID++ )
        {
          const auto vertexDoFIndex = neighboringVertex0 + elementAsIndices.at( localVertexID );
          macroCellStencilEntries[ vertexDoFIndex ] += weights[ localVertexID ];
        }
    }

    else {
      for ( uint_t localVertexID = 0; localVertexID < 4; localVertexID++ )
        {

          // find index into stencil
          const auto vertexDoFIndex = neighboringVertex0 + elementAsIndices.at( localVertexID );

          // obtain entry of local element matrix and add to stencil weight
          macroCellStencilEntries[ vertexDoFIndex ] += ufcGen.integrate( geometricCoordinates, centerEdge,
                                                                         { localVertexID, localVertexID } );
        }
    }

  }

  return macroCellStencilEntries;
}


/// \brief Calculates the edge to edge stencil weights for one orientation of the center edge and one of the leaves
///        from the stiffness matrices of neighboring elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-edge's location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
///
/// \param microEdgeIndex the logical index of the micro-edge in a macro-cell (can also be located on the macro-cell's boundary)
/// \param centerOrientation the orientation of the stencil center micro-edge
/// \param leafOrientation the orientation of the stencil leaves
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param ufcGen the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from index offsets to stencil weights
///
template< typename UFCOperator >
inline std::map< indexing::Index, real_t > calculateEdgeToEdgeStencilInMacroCell( const indexing::Index & microEdgeIndex,
                                                                                           const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                           const edgedof::EdgeDoFOrientation & leafOrientation,
                                                                                           const Cell & cell, const uint_t & level, const UFCOperator & ufcGen )
{
  typedef stencilDirection                     sd;
  std::map< indexing::Index, real_t > macroCellStencilEntries;

  const auto  offsetsToNeighborVertices = edgedof::calcNeighboringVertexDoFIndices( centerOrientation );
  const auto& neighboringVertex0        = offsetsToNeighborVertices.at( 0 );
  const auto  neighboringElementsAtVertex0 =
      P1Elements::P1Elements3D::getNeighboringElements( microEdgeIndex + neighboringVertex0, level );
  const auto secondDirectionInElements =
      vertexdof::stencilDirectionFromLogicalOffset( offsetsToNeighborVertices.at( 1 ) - offsetsToNeighborVertices.at( 0 ) );

  // 1. Going over all neighboring cells of the first neighboring micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions.
  for ( const auto& cellAtVertex : neighboringElementsAtVertex0 )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    // Exclude elements that do not contain the second vertex that defines the current edge.
    if ( std::find( cellAtVertex.begin(), cellAtVertex.end(), secondDirectionInElements ) == cellAtVertex.end() )
    {
        continue;
    }

    const std::array< indexing::Index, 4 > elementAsIndices = {
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[0] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[1] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[2] ),
        vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[3] ),
    };

    // Since there are no parallel edges in the same element, there is at most one edge with the specified orientation in the current element.
    // If there is none, we don't do anything with this element.
    if ( const auto optLeafEdge = edgeWithOrientationFromElement( elementAsIndices, leafOrientation ))
    {
        const auto leafEdge = optLeafEdge.value();

        // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex.
        //    The reference micro-index is the 'first' of the two indices that define the current edge.
        std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
        for ( uint_t localID = 0; localID < 4; localID++ )
        {
          logicalOffsetsFromCenter[localID] = microEdgeIndex + neighboringVertex0 + elementAsIndices.at( localID );
        }

        // 5. Adding contribution to stencil

        // obtain micro cell coordinates
        std::array< Point3D, 4 > geometricCoordinates;
        for ( uint_t localID = 0; localID < 4; localID++ )
        {
          geometricCoordinates[localID] =
              vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
        }

        // find index into stencil
      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( neighboringVertex0 + elementAsIndices.at( leafEdge.at(0) ),
                                                           neighboringVertex0 + elementAsIndices.at( leafEdge.at(1) ) );

      // find vertices of centerEdge
      auto centerEdge = edgeWithOrientationFromElement( elementAsIndices, centerOrientation ).value();

      // obtain entry of local element matrix and add to stencil weight
      macroCellStencilEntries[ edgeDoFIndex ] += ufcGen.integrate( geometricCoordinates, centerEdge, leafEdge );

    }
  }

  return macroCellStencilEntries;
}




}
}
}

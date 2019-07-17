
#pragma once

#include "core/Optional.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/p2functionspace/P2Form.hpp"
#include "tinyhhg_core/p2functionspace/generated_new/P2FenicsForm.hpp"

#include <set>

/// FENICS DoF ordering in 10x10 P2 stiffness matrix:
/// 0: vertex 0
/// 1: vertex 1
/// 2: vertex 2
/// 3: vertex 3
/// 4: edge (2, 3)
/// 5: edge (1, 3)
/// 6: edge (1, 2)
/// 7: edge (0, 3)
/// 8: edge (0, 2)
/// 9: edge (0, 1)

namespace hhg {
namespace P2Elements {
namespace P2Elements3D {

/// \brief Returns the (sorted) array indices of the two micro-vertices of the given element that span the edge with the specified orientation there are any.
inline walberla::optional< std::array< uint_t, 2 > > edgeWithOrientationFromElement( const std::array< indexing::IndexIncrement, 4 > & elementVertices,
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


/// \brief Calculates the local P2 stiffness matrix from the absolute indices of the element's vertices in a cell.
///
/// \param absoluteLogicalElementVertices the absolute logical indices of the element's vertices, their coordinates are given to the UFCOperator in the same order
/// \param cell the surrounding macro-cell
/// \param level the multigrid level
/// \param ufcGen the FENICS UFCOperator
/// \return the P2 stiffness matrix calculated from the given element and macro-cell
template< typename UFCOperator >
inline typename fenics::UFCTrait< UFCOperator >::LocalStiffnessMatrix_T calculateLocalP2StiffnessMatrix(
                                                  const std::array< indexing::Index, 4 > & absoluteLogicalElementVertices,
                                                  const Cell & cell, const uint_t & level, const UFCOperator & ufcGen )
{
  // Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
  std::array< Point3D, 4 > geometricCoordinates;
  for ( uint_t localID = 0; localID < 4; localID++ ) {
    geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, absoluteLogicalElementVertices[localID] );
  }

  std::array< Point3D, 4 > geometricOffsetsFromCenter;
  for ( uint_t localID = 0; localID < 4; localID++ ) {
    geometricOffsetsFromCenter[localID] = geometricCoordinates[localID] - geometricCoordinates[0];
  }

  // Computing the local stiffness matrix
  // To calculate the 10x10 stiffness matrix, we need the geometric offsets from the reference micro-vertex
  // from all micro-vertices in the neighbor cell (including the reference micro-vertex itself -> the first offset is always (0.0, 0.0, 0.0))

  // Flattening the offset array to be able to pass it to the fenics routines.
  double geometricOffsetsArray[12];
  for ( uint_t cellVertex = 0; cellVertex < 4; cellVertex++ ) {
    for ( uint_t coordinate = 0; coordinate < 3; coordinate++ ) {
      geometricOffsetsArray[cellVertex * 3 + coordinate] = geometricOffsetsFromCenter[cellVertex][coordinate];
    }
  }

  typename fenics::UFCTrait< UFCOperator >::LocalStiffnessMatrix_T localStiffnessMatrix;
  ufcGen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

  return localStiffnessMatrix;
}


/// \brief Calculates the edge to vertex stencil weights for one (leaf-)orientation from the stiffness matrices of neighboring elements at an index in a macro-cell.
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
inline std::map< indexing::IndexIncrement, real_t > calculateEdgeToVertexStencilInMacroCell( const indexing::Index & microVertexIndex,
                                                                                             const edgedof::EdgeDoFOrientation & leafOrientation,
                                                                                             const Cell & cell, const uint_t & level,
                                                                                             const UFCOperator & ufcGen )
{
  typedef stencilDirection sd;
  std::map< indexing::IndexIncrement, real_t > macroCellStencilEntries;

  const auto neighboringElements = P1Elements::P1Elements3D::getNeighboringElements( microVertexIndex, level );

  // 1. Going over all neighboring cells of a micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
  //    VERTEX_C represents the reference micro-vertex.
  for ( const auto & cellAtVertex : neighboringElements )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    const std::array< indexing::IndexIncrement, 4 > elementAsIndices = {
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

#ifdef OLD_SCHOOL
      const auto localStiffnessMatrix = calculateLocalP2StiffnessMatrix( logicalOffsetsFromCenter, cell, level, ufcGen );

      // 5. Adding contribution to stencil
      //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
      //    we always only need the first row of the local stiffness matrix.

      /// 4: edge (2, 3)
      /// 5: edge (1, 3)
      /// 6: edge (1, 2)
      /// 7: edge (0, 3)
      /// 8: edge (0, 2)
      /// 9: edge (0, 1)

      uint_t localEdgeIDInStiffnessMatrix = 10; // this should crash when accessing the matrix
           if ( edge.at(0) == 2 && edge.at(1) == 3 ) localEdgeIDInStiffnessMatrix = 4;
      else if ( edge.at(0) == 1 && edge.at(1) == 3 ) localEdgeIDInStiffnessMatrix = 5;
      else if ( edge.at(0) == 1 && edge.at(1) == 2 ) localEdgeIDInStiffnessMatrix = 6;
      else if ( edge.at(0) == 0 && edge.at(1) == 3 ) localEdgeIDInStiffnessMatrix = 7;
      else if ( edge.at(0) == 0 && edge.at(1) == 2 ) localEdgeIDInStiffnessMatrix = 8;
      else if ( edge.at(0) == 0 && edge.at(1) == 1 ) localEdgeIDInStiffnessMatrix = 9;
      else
      {
        WALBERLA_ASSERT( "Inconsistent element / orientation." );
      }

      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( elementAsIndices.at( edge.at(0) ), elementAsIndices.at( edge.at(1) ) );
      macroCellStencilEntries[ edgeDoFIndex ] += real_c( localStiffnessMatrix( 0, localEdgeIDInStiffnessMatrix ) );

#else
      // obtain micro cell coordinates
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ ) {
        geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      // find index into stencil
      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( elementAsIndices.at( edge.at(0) ), elementAsIndices.at( edge.at(1) ) );

      // obtain entry of local element matrix and add to stencil weight
      macroCellStencilEntries[ edgeDoFIndex ] += ufcGen.integrate( geometricCoordinates, { 0, 0 }, edge );

#endif

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
inline std::map< indexing::IndexIncrement, real_t > calculateVertexToEdgeStencilInMacroCell( const indexing::Index & microEdgeIndex,
                                                                                             const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                             const Cell & cell, const uint_t & level, const UFCOperator & ufcGen )
{
  typedef stencilDirection sd;
  std::map< indexing::IndexIncrement, real_t > macroCellStencilEntries;

  const auto offsetsToNeighborVertices = edgedof::calcNeighboringVertexDoFIndices( centerOrientation );
  const auto neighboringVertex0 = offsetsToNeighborVertices.at(0);
  const auto neighboringElementsAtVertex0 = P1Elements::P1Elements3D::getNeighboringElements( microEdgeIndex + neighboringVertex0, level );
  const auto secondDirectionInElements = vertexdof::stencilDirectionFromLogicalOffset( offsetsToNeighborVertices.at(1) - offsetsToNeighborVertices.at(0) );

  // 1. Going over all neighboring cells of the first neighboring micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions.
  for ( const auto & cellAtVertex : neighboringElementsAtVertex0 )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    // Exclude elements that do not contain the second vertex that defines the current edge.
    if ( std::find( cellAtVertex.begin(), cellAtVertex.end(), secondDirectionInElements ) == cellAtVertex.end() )
    {
      continue;
    }

    const std::array< indexing::IndexIncrement, 4 > elementAsIndices = {
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[0] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[1] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[2] ),
      vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[3] ),
    };

    // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex.
    //    The reference micro-index is the 'first' of the two indices that define the current edge.
    std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ ) {
      logicalOffsetsFromCenter[localID] = microEdgeIndex + neighboringVertex0 + elementAsIndices.at( localID );
    }

#ifdef OLD_SCHOOL
    const auto localStiffnessMatrix = calculateLocalP2StiffnessMatrix( logicalOffsetsFromCenter, cell, level, ufcGen );

    // 5. Adding contribution to stencil
    //    We now need to read from the correct row of the stiffness matrix.
    //    We know that the first vertex that defines the current edge is at position 0.
    //    Therefore (due to the ordering of in the FENICS matrix) the row is either at idx 7, 8 or 9.

    /// 4: edge (2, 3)
    /// 5: edge (1, 3)
    /// 6: edge (1, 2)
    /// 7: edge (0, 3)
    /// 8: edge (0, 2)
    /// 9: edge (0, 1)

    uint_t localEdgeIDInStiffnessMatrix = 10; // this should crash when accessing the matrix
    WALBERLA_ASSERT_NOT_IDENTICAL( cellAtVertex.at(0), secondDirectionInElements );
    if ( cellAtVertex.at(1) == secondDirectionInElements ) localEdgeIDInStiffnessMatrix = 9;
    else if ( cellAtVertex.at(2) == secondDirectionInElements ) localEdgeIDInStiffnessMatrix = 8;
    else if ( cellAtVertex.at(3) == secondDirectionInElements ) localEdgeIDInStiffnessMatrix = 7;
    else
    {
      WALBERLA_ASSERT( "Inconsistent element / orientation." );
    }

    for ( uint_t localVertexID = 0; localVertexID < 4; localVertexID++ )
    {
      const auto vertexDoFIndex = neighboringVertex0 + elementAsIndices.at( localVertexID );
      macroCellStencilEntries[ vertexDoFIndex ] += real_c( localStiffnessMatrix( localEdgeIDInStiffnessMatrix, localVertexID ) );
    }

#else
    // obtain micro cell coordinates
    std::array< Point3D, 4 > geometricCoordinates;
    for ( uint_t localID = 0; localID < 4; localID++ ) {
      geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
    }

    // find vertices of centerEdge
    auto centerEdge = edgeWithOrientationFromElement( elementAsIndices, centerOrientation ).value();

    for ( uint_t localVertexID = 0; localVertexID < 4; localVertexID++ )
    {

      // find index into stencil
      const auto vertexDoFIndex = neighboringVertex0 + elementAsIndices.at( localVertexID );

      // obtain entry of local element matrix and add to stencil weight
      macroCellStencilEntries[ vertexDoFIndex ] += ufcGen.integrate( geometricCoordinates, centerEdge,
                                                                     { localVertexID, localVertexID } );
    }

#endif
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
inline std::map< indexing::IndexIncrement, real_t > calculateEdgeToEdgeStencilInMacroCell( const indexing::Index & microEdgeIndex,
                                                                                           const edgedof::EdgeDoFOrientation & centerOrientation,
                                                                                           const edgedof::EdgeDoFOrientation & leafOrientation,
                                                                                           const Cell & cell, const uint_t & level, const UFCOperator & ufcGen )
{
  typedef stencilDirection sd;
  std::map< indexing::IndexIncrement, real_t > macroCellStencilEntries;

  const auto offsetsToNeighborVertices = edgedof::calcNeighboringVertexDoFIndices( centerOrientation );
  const auto neighboringVertex0 = offsetsToNeighborVertices.at(0);
  const auto neighboringElementsAtVertex0 = P1Elements::P1Elements3D::getNeighboringElements( microEdgeIndex + neighboringVertex0, level );
  const auto secondDirectionInElements = vertexdof::stencilDirectionFromLogicalOffset( offsetsToNeighborVertices.at(1) - offsetsToNeighborVertices.at(0) );

  // 1. Going over all neighboring cells of the first neighboring micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions.
  for ( const auto & cellAtVertex : neighboringElementsAtVertex0 )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    // Exclude elements that do not contain the second vertex that defines the current edge.
    if ( std::find( cellAtVertex.begin(), cellAtVertex.end(), secondDirectionInElements ) == cellAtVertex.end() )
    {
      continue;
    }

    const std::array< indexing::IndexIncrement, 4 > elementAsIndices = {
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
      for ( uint_t localID = 0; localID < 4; localID++ ) {
        logicalOffsetsFromCenter[localID] = microEdgeIndex + neighboringVertex0 + elementAsIndices.at( localID );
      }

#ifdef OLD_SCHOOL
      const auto localStiffnessMatrix = calculateLocalP2StiffnessMatrix( logicalOffsetsFromCenter, cell, level, ufcGen );

      // 5. Adding contribution to stencil
      //    We now need to read from the correct row of the stiffness matrix.
      //    We know that the first vertex that defines the current edge is at position 0.
      //    Therefore (due to the ordering of in the FENICS matrix) the row is either at idx 7, 8 or 9.

      /// 4: edge (2, 3)
      /// 5: edge (1, 3)
      /// 6: edge (1, 2)
      /// 7: edge (0, 3)
      /// 8: edge (0, 2)
      /// 9: edge (0, 1)

      uint_t centerEdgeIDInStiffnessMatrix = 10; // this should crash when accessing the matrix
      WALBERLA_ASSERT_NOT_IDENTICAL( cellAtVertex.at(0), secondDirectionInElements );
      if ( cellAtVertex.at(1) == secondDirectionInElements ) centerEdgeIDInStiffnessMatrix = 9;
      else if ( cellAtVertex.at(2) == secondDirectionInElements ) centerEdgeIDInStiffnessMatrix = 8;
      else if ( cellAtVertex.at(3) == secondDirectionInElements ) centerEdgeIDInStiffnessMatrix = 7;
      else
      {
        WALBERLA_ASSERT( "Inconsistent element / orientation." );
      }

      uint_t leafEdgeIDInStiffnessMatrix = 10; // this should crash when accessing the matrix
      if ( leafEdge.at(0) == 2 && leafEdge.at(1) == 3 ) leafEdgeIDInStiffnessMatrix = 4;
      else if ( leafEdge.at(0) == 1 && leafEdge.at(1) == 3 ) leafEdgeIDInStiffnessMatrix = 5;
      else if ( leafEdge.at(0) == 1 && leafEdge.at(1) == 2 ) leafEdgeIDInStiffnessMatrix = 6;
      else if ( leafEdge.at(0) == 0 && leafEdge.at(1) == 3 ) leafEdgeIDInStiffnessMatrix = 7;
      else if ( leafEdge.at(0) == 0 && leafEdge.at(1) == 2 ) leafEdgeIDInStiffnessMatrix = 8;
      else if ( leafEdge.at(0) == 0 && leafEdge.at(1) == 1 ) leafEdgeIDInStiffnessMatrix = 9;
      else
      {
        WALBERLA_ASSERT( "Inconsistent element / orientation." );
      }

      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( neighboringVertex0 + elementAsIndices.at( leafEdge.at(0) ), neighboringVertex0 + elementAsIndices.at( leafEdge.at(1) ) );
      macroCellStencilEntries[ edgeDoFIndex ] += real_c( localStiffnessMatrix( centerEdgeIDInStiffnessMatrix, leafEdgeIDInStiffnessMatrix ) );

#else

      // obtain micro cell coordinates
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ ) {
        geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      // find index into stencil
      const auto edgeDoFIndex = edgedof::calcEdgeDoFIndex( neighboringVertex0 + elementAsIndices.at( leafEdge.at(0) ),
                                                           neighboringVertex0 + elementAsIndices.at( leafEdge.at(1) ) );

      // find vertices of centerEdge
      auto centerEdge = edgeWithOrientationFromElement( elementAsIndices, centerOrientation ).value();

      // obtain entry of local element matrix and add to stencil weight
      macroCellStencilEntries[ edgeDoFIndex ] += ufcGen.integrate( geometricCoordinates, centerEdge, leafEdge );

#endif

    }
  }

  return macroCellStencilEntries;
}




}
}
}

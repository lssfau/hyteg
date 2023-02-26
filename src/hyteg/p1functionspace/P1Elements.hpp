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

#include <unordered_map>

#include "hyteg/StencilDirections.hpp"
#include "hyteg/fenics/fenics.hpp"
#include "hyteg/fenics/ufc_traits.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

namespace P1Elements {
namespace P1Elements2D {

// Fenics P1 DoF ordering
// 2       1---0
// |\       \  |
// | \       \ |
// |  \       \|
// 0---1       2

const uint_t ElementSize = 3;

typedef stencilDirection                  SD;
typedef std::array< SD, ElementSize >     P1Element;
typedef std::array< uint_t, ElementSize > DoFMap;
typedef std::array< uint_t, ElementSize > StencilMap;

const P1Element elementSW = { { SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S } };
const P1Element elementS  = { { SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE } };
const P1Element elementSE = { { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E } };
const P1Element elementNE = { { SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N } };
const P1Element elementN  = { { SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_NW } };
const P1Element elementNW = { { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_W } };

inline StencilMap convertStencilDirectionsToIndices( const P1Element& element )
{
   return { { vertexdof::stencilIndexFromVertex( element[0] ),
              vertexdof::stencilIndexFromVertex( element[1] ),
              vertexdof::stencilIndexFromVertex( element[2] ) } };
}

template < typename StencilMemory >
inline void assembleP1LocalStencil( const StencilMap& stencilMap,
                                    const DoFMap&     dofMap,
                                    const Matrix3r&   localMatrix,
                                    StencilMemory&    stencil,
                                    double            coeffWeight = 1.0 )
{
   for ( uint_t j = 0; j < 3; ++j )
   {
      stencil[stencilMap[j]] += coeffWeight * localMatrix( dofMap[0], dofMap[j] );
   }
}

} // namespace P1Elements2D

namespace P1Elements3D {

typedef stencilDirection sd;

const std::array< std::array< stencilDirection, 4 >, 4 > whiteUpCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN }, // below
    { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },  // top front
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },  // top back west
    { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },   // top back east
} };

const std::array< std::array< stencilDirection, 4 >, 4 > whiteDownCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_S },   // below front west
    { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_BE },  // below front east
    { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_NW, sd::VERTEX_BN },  // below back
    { sd::VERTEX_C, sd::VERTEX_TS, sd::VERTEX_TC, sd::VERTEX_TW }, // top
} };

const std::array< std::array< stencilDirection, 4 >, 4 > blueUpCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BN, sd::VERTEX_BNW }, // below
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_S, sd::VERTEX_TS },    // top front west
    { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_TSE },  // top front east
    { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC },   // top back
} };

const std::array< std::array< stencilDirection, 4 >, 4 > blueDownCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE },   // below front
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW },  // below back west
    { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N },    // below back east
    { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE }, // top
} };

const std::array< std::array< stencilDirection, 4 >, 4 > greenUpCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW },  // below west
    { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN },   // below east
    { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW },  // top back
    { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE }, // top front
} };

const std::array< std::array< stencilDirection, 4 >, 4 > greenDownCellsAtInnerVertex = { {
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE },  // below front
    { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW }, // below back
    { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC },  // top east
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW },   // top west
} };

const std::array< std::array< stencilDirection, 4 >, 24 > allCellsAtInnerVertex = {
    { whiteUpCellsAtInnerVertex[0],   whiteUpCellsAtInnerVertex[1],   whiteUpCellsAtInnerVertex[2],
      whiteUpCellsAtInnerVertex[3],   whiteDownCellsAtInnerVertex[0], whiteDownCellsAtInnerVertex[1],
      whiteDownCellsAtInnerVertex[2], whiteDownCellsAtInnerVertex[3], blueUpCellsAtInnerVertex[0],
      blueUpCellsAtInnerVertex[1],    blueUpCellsAtInnerVertex[2],    blueUpCellsAtInnerVertex[3],
      blueDownCellsAtInnerVertex[0],  blueDownCellsAtInnerVertex[1],  blueDownCellsAtInnerVertex[2],
      blueDownCellsAtInnerVertex[3],  greenUpCellsAtInnerVertex[0],   greenUpCellsAtInnerVertex[1],
      greenUpCellsAtInnerVertex[2],   greenUpCellsAtInnerVertex[3],   greenDownCellsAtInnerVertex[0],
      greenDownCellsAtInnerVertex[1], greenDownCellsAtInnerVertex[2], greenDownCellsAtInnerVertex[3] } };

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace0 = {
    { // no cells with bottom direction
      { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
      { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_TS, sd::VERTEX_TC, sd::VERTEX_TW },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_S, sd::VERTEX_TS },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW },
      { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW } } };

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace1 = {
    { // no cells with south direction
      { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
      { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_NW, sd::VERTEX_BN },
      { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BN, sd::VERTEX_BNW },
      { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N },
      { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN },
      { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW },
      { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW } } };

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace2 = {
    { // no cells with west direction
      { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
      { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
      { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_BE },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N },
      { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN },
      { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE },
      { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE },
      { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC } } };

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace3 = { {
    // no cells in {N, E, TC, TSE}
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
    { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_S },
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BN, sd::VERTEX_BNW },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_S, sd::VERTEX_TS },
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW },
    { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE },
    { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW },
    { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW },
} };

const std::array< std::array< std::array< stencilDirection, 4 >, 12 >, 4 > allCellsAtFace = {
    { allCellsAtFace0, allCellsAtFace1, allCellsAtFace2, allCellsAtFace3 } };

/// \brief Returns the neighboring elements (== micro-cells) of a given micro-vertex index
///
/// \param microVertexIndex the micro-vertex index
/// \param level the multigrid level
/// \return a vector filled with 4-tuples (std::array) of stencil directions from that micro-vertex index that construct a micro-cell each (therefore always contains VERTEX_C)
///
inline std::vector< std::array< stencilDirection, 4 > > getNeighboringElements( const indexing::Index& microVertexIndex,
                                                                                const uint_t&          level )
{
   typedef std::vector< std::array< stencilDirection, 4 > > returnType;

   WALBERLA_ASSERT_LESS( microVertexIndex.x() + microVertexIndex.y() + microVertexIndex.z(),
                         levelinfo::num_microvertices_per_edge( level ) );

   const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
   const auto onCellEdges    = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
   const auto onCellFaces    = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );

   if ( level == 0 )
   {
      WALBERLA_ASSERT( !onCellVertices.empty() );
   }

   if ( onCellVertices.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( onCellVertices.size(), 1 );
      WALBERLA_ASSERT_EQUAL( onCellEdges.size(), 3 );
      WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 3 );
      const auto localVertexID   = *onCellVertices.begin();
      const auto singleMicroCell = [localVertexID] {
         switch ( localVertexID )
         {
         case 0:
            return std::array< stencilDirection, 4 >( { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC } );
         case 1:
            return std::array< stencilDirection, 4 >( { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW } );
         case 2:
            return std::array< stencilDirection, 4 >( { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS } );
         default:
            return std::array< stencilDirection, 4 >( { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN } );
         }
      }();
      return returnType( { singleMicroCell } );
   }
   else if ( onCellEdges.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( onCellEdges.size(), 1 );
      const auto localEdgeID = *onCellEdges.begin();
      WALBERLA_ASSERT_GREATER_EQUAL( localEdgeID, 0 );
      WALBERLA_ASSERT_LESS_EQUAL( localEdgeID, 6 );
      WALBERLA_UNUSED( localEdgeID );

      WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 2 );
      const std::vector< uint_t > onCellFacesVector( onCellFaces.begin(), onCellFaces.end() );

      auto cellsAtFace0 = allCellsAtFace[onCellFacesVector[0]];
      auto cellsAtFace1 = allCellsAtFace[onCellFacesVector[1]];
      std::sort( cellsAtFace0.begin(), cellsAtFace0.end() );
      std::sort( cellsAtFace1.begin(), cellsAtFace1.end() );

      returnType allCellsAtEdge;
      std::set_intersection( cellsAtFace0.begin(),
                             cellsAtFace0.end(),
                             cellsAtFace1.begin(),
                             cellsAtFace1.end(),
                             std::back_inserter( allCellsAtEdge ) );
      return allCellsAtEdge;
   }
   else if ( onCellFaces.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 1 );
      const auto localFaceID = *onCellFaces.begin();
      WALBERLA_ASSERT_GREATER_EQUAL( localFaceID, 0 );
      WALBERLA_ASSERT_LESS_EQUAL( localFaceID, 3 );
      return returnType( allCellsAtFace[localFaceID].begin(), allCellsAtFace[localFaceID].end() );
   }
   else
   {
      return returnType( allCellsAtInnerVertex.begin(), allCellsAtInnerVertex.end() );
   }
}

/// \brief Calculates the stencil weights from the stiffness matrices of neighboring elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-vertex' location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
///
/// \param microVertexIndex the logical index of the micro-vertex in a macro-cell (can also be located on the macro-cell's boundary)
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param ufcGen the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from stencil directions to stencil weights
///
template < typename UFCOperator >
inline std::map< stencilDirection, real_t > calculateStencilInMacroCell( const indexing::Index& microVertexIndex,
                                                                         const Cell&            cell,
                                                                         const uint_t&          level,
                                                                         const UFCOperator&     ufcGen )
{
   std::map< stencilDirection, real_t > macroCellStencilEntries;

   const auto neighboringElements = getNeighboringElements( microVertexIndex, level );

   // 1. Going over all neighboring cells of a micro-vertex
   //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
   //    VERTEX_C represents the reference micro-vertex.
   for ( const auto& cellAtVertex : neighboringElements )
   {
      WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

      // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
      std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         logicalOffsetsFromCenter[localID] =
             microVertexIndex + indexing::Index( vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[localID] ) );
      }

      // 3. Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         geometricCoordinates[localID] =
             vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      std::array< Point3D, 4 > geometricOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         geometricOffsetsFromCenter[localID] = geometricCoordinates[localID] - geometricCoordinates[0];
      }

      // 4. Computing the local stiffness matrix
      //    To calculate the 4x4 stiffness matrix, we need the geometric offsets from the reference micro-vertex
      //    from all micro-vertices in the neighbor cell (including the reference micro-vertex itself -> the first offset is always (0.0, 0.0, 0.0))

      // Flattening the offset array to be able to pass it to the fenics routines.
      double geometricOffsetsArray[12];
      for ( uint_t cellVertex = 0; cellVertex < 4; cellVertex++ )
      {
         for ( uint_t coordinate = 0; coordinate < 3; coordinate++ )
         {
            geometricOffsetsArray[cellVertex * 3 + coordinate] = geometricOffsetsFromCenter[cellVertex][coordinate];
         }
      }

      typename fenics::UFCTrait< UFCOperator >::LocalStiffnessMatrix_T localStiffnessMatrix;
      ufcGen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

      // 5. Adding contribution to stencil
      //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
      //    we always only need the first row of the local stiffness matrix.
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         const stencilDirection stencilDir = cellAtVertex[localID];
         if ( macroCellStencilEntries.count( stencilDir ) == 0 )
         {
            macroCellStencilEntries[stencilDir] = real_c( 0 );
         }
         macroCellStencilEntries[stencilDir] += real_c( localStiffnessMatrix( 0, localID ) );
      }
   }
   return macroCellStencilEntries;
}

/// \brief Calculates the stencil weights from the stiffness matrices of neighboring elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-vertex' location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
/// USE WITH CAUTION: ONLY WORKS FOR CONSTANT STENCILS!
/// \param microVertexIndex the logical index of the micro-vertex in a macro-cell (can also be located on the macro-cell's boundary)
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from stencil directions to stencil weights
///
template < class P1Form >
inline std::map< stencilDirection, real_t > calculateStencilInMacroCellForm( const indexing::Index& microVertexIndex,
                                                                             const Cell&            cell,
                                                                             const uint_t&          level,
                                                                             const P1Form&          form )
{
   std::map< stencilDirection, real_t > macroCellStencilEntries;

   const auto neighboringElements = getNeighboringElements( microVertexIndex, level );

   // 1. Going over all neighboring cells of a micro-vertex
   //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
   //    VERTEX_C represents the reference micro-vertex.
   for ( const auto& cellAtVertex : neighboringElements )
   {
      WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

      // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
      std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         logicalOffsetsFromCenter[localID] =
             microVertexIndex + indexing::Index( vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[localID] ) );
      }

      // 3. Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         geometricCoordinates[localID] =
             vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      std::array< Point3D, 4 > geometricOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         geometricOffsetsFromCenter[localID] = geometricCoordinates[localID] - geometricCoordinates[0];
      }

      // 4. Computing the local stiffness matrix
      //    To calculate the 4x4 stiffness matrix, we need the geometric offsets from the reference micro-vertex
      //    from all micro-vertices in the neighbor cell (including the reference micro-vertex itself -> the first offset is always (0.0, 0.0, 0.0))
      //!!! This shift results in wrong stencils unless the stencil is constant !!!
      Point4D localStiffnessMatrixRow;
      form.integrate( geometricOffsetsFromCenter, localStiffnessMatrixRow );

      // 5. Adding contribution to stencil
      //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
      //    we only need to get the result of the form integrator which gives us the first row of the local stiffness matrix
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         const stencilDirection stencilDir = cellAtVertex[localID];
         if ( macroCellStencilEntries.count( stencilDir ) == 0 )
         {
            macroCellStencilEntries[stencilDir] = real_c( 0 );
         }
         macroCellStencilEntries[stencilDir] += real_c( localStiffnessMatrixRow[localID] );
      }
   }
   return macroCellStencilEntries;
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
/// \brief Calculates the stencil weights from the stiffness matrices of neighboring elements at an index in a macro-cell.
///
/// Also works for indices on the boundary of a macro-cell. In this case the stencil map simply contains less elements.
/// It automatically computes / selects the neighboring elements depending on the micro-vertex' location.
/// Note that only the weights for the stencil that lie in the specified macro-cell are returned.
///
/// \param microVertexIndex the logical index of the micro-vertex in a macro-cell (can also be located on the macro-cell's boundary)
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a (variable sized) map from stencil directions to stencil weights
///
template < class P1Form >
inline std::map< stencilDirection, real_t > calculateStencilInMacroCellForm_new( const indexing::Index& microVertexIndex,
                                                                                 const Cell&            cell,
                                                                                 const uint_t&          level,
                                                                                 P1Form&                form )
{
   form.setGeometryMap( cell.getGeometryMap() );

   std::map< stencilDirection, real_t > macroCellStencilEntries;

   const auto neighboringElements = getNeighboringElements( microVertexIndex, level );

   // 1. Going over all neighboring cells of a micro-vertex
   //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
   //    VERTEX_C represents the reference micro-vertex.
   for ( const auto& cellAtVertex : neighboringElements )
   {
      WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

      // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
      std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         logicalOffsetsFromCenter[localID] = microVertexIndex + vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[localID] );
      }

      // 3. Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
      std::array< Point3D, 4 > geometricCoordinates;
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         geometricCoordinates[localID] =
             vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
      }

      // 4. Computing the local stiffness matrix
      Matrixr< 1, 4 > localStiffnessMatrixRow;
      form.integrateRow( 0, geometricCoordinates, localStiffnessMatrixRow );
      // std::cout << " " << localStiffnessMatrixRow(0,0);

      // 5. Adding contribution to stencil
      //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
      //    we only need to get the result of the form integrator which gives us the first row of the local stiffness matrix
      for ( uint_t localID = 0; localID < 4; localID++ )
      {
         const stencilDirection stencilDir = cellAtVertex[localID];
         if ( macroCellStencilEntries.count( stencilDir ) == 0 )
         {
            macroCellStencilEntries[stencilDir] = real_c( 0 );
         }
         macroCellStencilEntries[stencilDir] += real_c( localStiffnessMatrixRow( 0, localID ) );
      }
   }
   return macroCellStencilEntries;
}

/// \brief Assembles the local P1 operator stencil on a macro-vertex
///
/// \param storage the governing \ref PrimitiveStorage
/// \param vertex the macro-vertex
/// \param microVertexIndex the micro-vertex index on the macro-vertex (currently this must be (0, 0, 0))
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a vector containing the stencil weights for the micro-vertex on that macro-vertex,
///         stencil[0] is the center weight, stencil[neighborID + 1] is the weight for the neighbor with neighborID
///
template < class P1Form >
inline std::vector< real_t > assembleP1LocalStencil( const std::shared_ptr< PrimitiveStorage >& storage,
                                                     const Vertex&                              vertex,
                                                     const indexing::Index&                     microVertexIndex,
                                                     const uint_t&                              level,
                                                     const P1Form&                              form )
{
   WALBERLA_CHECK_EQUAL(
       microVertexIndex, indexing::Index( 0, 0, 0 ), "[P1 vertex stencil assembly] micro-vertex index must be (0, 0, 0)" );

   const uint_t          stencilSize = vertexDoFMacroVertexStencilMemorySize( level, vertex );
   std::vector< real_t > stencil( stencilSize, real_c( 0 ) );

   for ( const auto& macroCellID : vertex.neighborCells() )
   {
      const auto macroCell = storage->getCell( macroCellID );

      // 1. translate coordinate to macro-cell
      const uint_t localVertexID = macroCell->getLocalVertexID( vertex.getID() );
      WALBERLA_ASSERT_LESS_EQUAL( localVertexID, 3 );
      const indexing::Index indexInMacroCell = [localVertexID, level] {
         switch ( localVertexID )
         {
         case 0:
            return indexing::Index( 0, 0, 0 );
         case 1:
            return indexing::Index( levelinfo::num_microvertices_per_edge( level ) - 1, 0, 0 );
         case 2:
            return indexing::Index( 0, levelinfo::num_microvertices_per_edge( level ) - 1, 0 );
         default:
            return indexing::Index( 0, 0, levelinfo::num_microvertices_per_edge( level ) - 1 );
         }
      }();

      WALBERLA_DEBUG_SECTION()
      {
         const auto debugLocalVertices = vertexdof::macrocell::isOnCellVertex( indexInMacroCell, level );
         const auto debugLocalEdges    = vertexdof::macrocell::isOnCellEdge( indexInMacroCell, level );
         const auto debugLocalFaces    = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         WALBERLA_ASSERT_EQUAL( debugLocalVertices.size(), 1 );
         WALBERLA_ASSERT_EQUAL( debugLocalEdges.size(), 3 );
         WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 3 );
         WALBERLA_ASSERT_EQUAL( *debugLocalVertices.begin(), localVertexID );
      }

      // 2. calculate stiffness matrix for each micro-cell and store contributions
      const auto cellLocalStencilWeights = calculateStencilInMacroCellForm( indexInMacroCell, *macroCell, level, form );

      // 3. translate coordinates / stencil directions back to vertex-local coordinate system
      for ( const auto it : cellLocalStencilWeights )
      {
         const auto            cellLocalDir        = it.first;
         const auto            stencilWeight       = it.second;
         const indexing::Index cellLocalIndexInDir = indexInMacroCell + vertexdof::logicalIndexOffsetFromVertex( cellLocalDir );
         const auto            onLocalEdgesDir     = vertexdof::macrocell::isOnCellEdge( cellLocalIndexInDir, level );
         const auto            onLocalVerticesDir  = vertexdof::macrocell::isOnCellVertex( cellLocalIndexInDir, level );
         if ( onLocalEdgesDir.size() == 1 )
         {
            const auto cellLocalEdgeID   = *onLocalEdgesDir.begin();
            const auto edgePrimitiveID   = macroCell->neighborEdges()[cellLocalEdgeID];
            const auto vertexLocalEdgeID = vertex.edge_index( edgePrimitiveID );
            stencil[vertexLocalEdgeID + 1] += stencilWeight;
         }
         else if ( onLocalVerticesDir.size() == 1 && level == 0 && cellLocalDir != sd::VERTEX_C )
         {
            const auto cellLocalVertexIDOfLeaf = *onLocalVerticesDir.begin();
            const auto cellLocalEdgeID =
                indexing::getCellLocalEdgeIDFromCellLocalVertexIDs( localVertexID, cellLocalVertexIDOfLeaf );
            const auto edgePrimitiveID   = macroCell->neighborEdges()[cellLocalEdgeID];
            const auto vertexLocalEdgeID = vertex.edge_index( edgePrimitiveID );
            stencil[vertexLocalEdgeID + 1] += stencilWeight;
         }
         else
         {
            WALBERLA_ASSERT_EQUAL( onLocalEdgesDir.size(), 3 );
            WALBERLA_ASSERT_EQUAL( cellLocalDir, sd::VERTEX_C );
            stencil[0] += stencilWeight;
         }
      }
   }

   return stencil;
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
/// \brief Assembles the local P1 operator stencil on a macro-vertex
///
/// \param storage the governing \ref PrimitiveStorage
/// \param vertex the macro-vertex
/// \param microVertexIndex the micro-vertex index on the macro-vertex (currently this must be (0, 0, 0))
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a vector containing the stencil weights for the micro-vertex on that macro-vertex,
///         stencil[0] is the center weight, stencil[neighborID + 1] is the weight for the neighbor with neighborID
///
template < class P1Form >
inline std::vector< real_t > assembleP1LocalStencil_new( const std::shared_ptr< PrimitiveStorage >& storage,
                                                         const Vertex&                              vertex,
                                                         const indexing::Index&                     microVertexIndex,
                                                         const uint_t&                              level,
                                                         P1Form&                                    form )
{
   WALBERLA_CHECK_EQUAL(
       microVertexIndex, indexing::Index( 0, 0, 0 ), "[P1 vertex stencil assembly] micro-vertex index must be (0, 0, 0)" );

   const uint_t          stencilSize = vertexDoFMacroVertexStencilMemorySize( level, vertex );
   std::vector< real_t > stencil( stencilSize, real_c( 0 ) );

   for ( const auto& macroCellID : vertex.neighborCells() )
   {
      const auto macroCell = storage->getCell( macroCellID );

      // 1. translate coordinate to macro-cell
      const uint_t localVertexID = macroCell->getLocalVertexID( vertex.getID() );
      WALBERLA_ASSERT_LESS_EQUAL( localVertexID, 3 );
      const indexing::Index indexInMacroCell = [localVertexID, level] {
         switch ( localVertexID )
         {
         case 0:
            return indexing::Index( 0, 0, 0 );
         case 1:
            return indexing::Index( levelinfo::num_microvertices_per_edge( level ) - 1, 0, 0 );
         case 2:
            return indexing::Index( 0, levelinfo::num_microvertices_per_edge( level ) - 1, 0 );
         default:
            return indexing::Index( 0, 0, levelinfo::num_microvertices_per_edge( level ) - 1 );
         }
      }();

      WALBERLA_DEBUG_SECTION()
      {
         const auto debugLocalVertices = vertexdof::macrocell::isOnCellVertex( indexInMacroCell, level );
         const auto debugLocalEdges    = vertexdof::macrocell::isOnCellEdge( indexInMacroCell, level );
         const auto debugLocalFaces    = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         WALBERLA_ASSERT_EQUAL( debugLocalVertices.size(), 1 );
         WALBERLA_ASSERT_EQUAL( debugLocalEdges.size(), 3 );
         WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 3 );
         WALBERLA_ASSERT_EQUAL( *debugLocalVertices.begin(), localVertexID );
      }

      // 2. calculate stiffness matrix for each micro-cell and store contributions
      const auto cellLocalStencilWeights = calculateStencilInMacroCellForm_new( indexInMacroCell, *macroCell, level, form );

      // 3. translate coordinates / stencil directions back to vertex-local coordinate system
      for ( const auto it : cellLocalStencilWeights )
      {
         const auto                            cellLocalDir  = it.first;
         const auto                            stencilWeight = it.second;
         const hyteg::indexing::Index cellLocalIndexInDir =
             indexInMacroCell + vertexdof::logicalIndexOffsetFromVertex( cellLocalDir );
         const auto onLocalEdgesDir    = vertexdof::macrocell::isOnCellEdge( cellLocalIndexInDir, level );
         const auto onLocalVerticesDir = vertexdof::macrocell::isOnCellVertex( cellLocalIndexInDir, level );
         if ( onLocalEdgesDir.size() == 1 )
         {
            const auto cellLocalEdgeID   = *onLocalEdgesDir.begin();
            const auto edgePrimitiveID   = macroCell->neighborEdges()[cellLocalEdgeID];
            const auto vertexLocalEdgeID = vertex.edge_index( edgePrimitiveID );
            stencil[vertexLocalEdgeID + 1] += stencilWeight;
         }
         else if ( onLocalVerticesDir.size() == 1 && level == 0 && cellLocalDir != sd::VERTEX_C )
         {
            const auto cellLocalVertexIDOfLeaf = *onLocalVerticesDir.begin();
            const auto cellLocalEdgeID =
                indexing::getCellLocalEdgeIDFromCellLocalVertexIDs( localVertexID, cellLocalVertexIDOfLeaf );
            const auto edgePrimitiveID   = macroCell->neighborEdges()[cellLocalEdgeID];
            const auto vertexLocalEdgeID = vertex.edge_index( edgePrimitiveID );
            stencil[vertexLocalEdgeID + 1] += stencilWeight;
         }
         else
         {
            WALBERLA_ASSERT_EQUAL( onLocalEdgesDir.size(), 3 );
            WALBERLA_ASSERT_EQUAL( cellLocalDir, sd::VERTEX_C );
            stencil[0] += stencilWeight;
         }
      }
   }

   return stencil;
}

/// \brief Assembles the local P1 operator stencil on a macro-edge
///
/// \param storage the governing \ref PrimitiveStorage
/// \param edge the macro-edge
/// \param microVertexIndex the micro-vertex index on the macro-edge (the y and z coordinate must be 0)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a vector containing the stencil weights for the micro-vertex on that macro-edge,
///         weights are sorted according to the vertexdof-macro-edge stencil index function
///
template < class P1Form >
inline std::vector< real_t > assembleP1LocalStencil( const std::shared_ptr< PrimitiveStorage >& storage,
                                                     const Edge&                                edge,
                                                     const indexing::Index&                     microVertexIndex,
                                                     const uint_t&                              level,
                                                     const P1Form&                              form )
{
   // check if index lies in the edges's interior
   WALBERLA_CHECK_EQUAL( microVertexIndex.y(), 0, "[P1 edge stencil assembly] y-coordinate on edge must be zero" );
   WALBERLA_CHECK_EQUAL( microVertexIndex.z(), 0, "[P1 edge stencil assembly] z-coordinate on edge must be zero" );
   WALBERLA_CHECK_GREATER( microVertexIndex.x(), 0 );
   WALBERLA_CHECK_LESS( microVertexIndex.x(), levelinfo::num_microvertices_per_edge( level ) );

   const uint_t          stencilSize = vertexDoFMacroEdgeStencilMemorySize( level, edge );
   std::vector< real_t > stencil( stencilSize, real_c( 0 ) );

   if ( level == 0 )
      return stencil;

   for ( const auto& macroCellID : edge.neighborCells() )
   {
      const auto macroCell = storage->getCell( macroCellID );

      // 1. translate coordinate to macro-cell

      // find out the local ID of the edge in the cell
      const uint_t localEdgeID = macroCell->getLocalEdgeID( edge.getID() );

      // Find out the coordinate system basis of the index on the macro-cell.
      WALBERLA_ASSERT_EQUAL( macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].size(), 2 );
      const uint_t basisCenter     = macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 0 );
      const uint_t basisXDirection = macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 1 );
      // find out the missing Z direction
      const std::set< uint_t > allDirections      = { 0, 1, 2, 3 };
      const std::set< uint_t > allDirectionsButYZ = { basisCenter, basisXDirection };
      std::vector< uint_t >    missingDirections;
      std::set_difference( allDirections.begin(),
                           allDirections.end(),
                           allDirectionsButYZ.begin(),
                           allDirectionsButYZ.end(),
                           std::inserter( missingDirections, missingDirections.begin() ) );
      WALBERLA_ASSERT_EQUAL( missingDirections.size(), 2 );
      const uint_t                  basisYDirection  = missingDirections[0];
      const uint_t                  basisZDirection  = missingDirections[1];
      const std::array< uint_t, 4 > indexingBasis    = { basisCenter, basisXDirection, basisYDirection, basisZDirection };
      const auto                    indexInMacroCell = indexing::basisConversion(
          microVertexIndex, indexingBasis, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );

      WALBERLA_DEBUG_SECTION()
      {
         const auto debugLocalEdges = vertexdof::macrocell::isOnCellEdge( indexInMacroCell, level );
         const auto debugLocalFaces = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         WALBERLA_ASSERT_EQUAL( debugLocalEdges.size(), 1 );
         WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 2 );
         WALBERLA_ASSERT_EQUAL( *debugLocalEdges.begin(), localEdgeID );
      }

      // 2. calculate stiffness matrix for each micro-cell and store contributions
      const auto cellLocalStencilWeights = calculateStencilInMacroCellForm( indexInMacroCell, *macroCell, level, form );

      // 3. translate coordinates / stencil directions back to edge-local coordinate system
      for ( const auto it : cellLocalStencilWeights )
      {
         const auto cellLocalDir  = it.first;
         const auto stencilWeight = it.second;

         const indexing::Index cellLocalIndexInDir = indexInMacroCell + vertexdof::logicalIndexOffsetFromVertex( cellLocalDir );
         const auto            onLocalFacesCenter  = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         const auto            onLocalFacesDir     = vertexdof::macrocell::isOnCellFace( cellLocalIndexInDir, level );
         std::vector< uint_t > intersectingFaces;
         std::set_intersection( onLocalFacesCenter.begin(),
                                onLocalFacesCenter.end(),
                                onLocalFacesDir.begin(),
                                onLocalFacesDir.end(),
                                std::back_inserter( intersectingFaces ) );

         if ( intersectingFaces.size() >= 2 )
         {
            // on edge
            const auto edgeLocalIndexInDir = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, indexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            WALBERLA_ASSERT_EQUAL( edgeLocalIndexInDir.y(), 0 );
            WALBERLA_ASSERT_EQUAL( edgeLocalIndexInDir.z(), 0 );
            const int dirDIfference = static_cast< int >( edgeLocalIndexInDir.x() - microVertexIndex.x() );
            WALBERLA_ASSERT_GREATER_EQUAL( dirDIfference, -1 );
            WALBERLA_ASSERT_LESS_EQUAL( dirDIfference, 1 );
            const stencilDirection dirOnEdge =
                dirDIfference == 0 ? sd::VERTEX_C : ( dirDIfference == 1 ? sd::VERTEX_E : sd::VERTEX_W );
            stencil[vertexdof::macroedge::stencilIndexOnEdge( dirOnEdge )] += stencilWeight;
         }
         else if ( intersectingFaces.size() == 1 )
         {
            // on neighbor face
            const auto localFaceIDInCell = *intersectingFaces.begin();
            const auto facePrimitiveID   = macroCell->neighborFaces()[localFaceIDInCell];
            // To get the correct indexing basis, we check which one results in a zero entry in the z coordinate.
            const auto                    firstTestIndexingBasis  = indexingBasis;
            const std::array< uint_t, 4 > secondTestIndexingBasis = {
                indexingBasis[0], indexingBasis[1], indexingBasis[3], indexingBasis[2] };
            const auto edgeLocalIndexInDirFirst = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, firstTestIndexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            const auto edgeLocalIndexInDirSecond = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, secondTestIndexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            WALBERLA_ASSERT_UNEQUAL( edgeLocalIndexInDirFirst.z(), edgeLocalIndexInDirSecond.z() );
            WALBERLA_ASSERT( edgeLocalIndexInDirFirst.z() == 0 || edgeLocalIndexInDirSecond.z() == 0 );
            const hyteg::indexing::Index faceLocalIndexInDir =
                edgeLocalIndexInDirFirst.z() == 0 ? edgeLocalIndexInDirFirst : edgeLocalIndexInDirSecond;
            WALBERLA_ASSERT_EQUAL( faceLocalIndexInDir.y(), 1 );
            stencilDirection faceLocalStencilDirection;

            const auto xOffset = static_cast< int >( faceLocalIndexInDir.x() ) - static_cast< int >( microVertexIndex.x() );
            if ( xOffset == 0 )
            {
               faceLocalStencilDirection = sd::VERTEX_E;
            }
            else if ( xOffset == -1 )
            {
               faceLocalStencilDirection = sd::VERTEX_W;
            }
            else
            {
               WALBERLA_ABORT( "[P1Elements][Edge] Invalid offsets" << xOffset );
            }
            stencil[vertexdof::macroedge::stencilIndexOnNeighborFace( faceLocalStencilDirection,
                                                                      edge.face_index( facePrimitiveID ) )] += stencilWeight;
         }
         else if ( intersectingFaces.size() == 0 )
         {
            // in macro-cell
            stencil[vertexdof::macroedge::stencilIndexOnNeighborCell( edge.cell_index( macroCellID ),
                                                                      edge.getNumNeighborFaces() )] += stencilWeight;
         }
      }
   }
   return stencil;
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
/// \brief Assembles the local P1 operator stencil on a macro-edge
///
/// \param storage the governing \ref PrimitiveStorage
/// \param edge the macro-edge
/// \param microVertexIndex the micro-vertex index on the macro-edge (the y and z coordinate must be 0)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a vector containing the stencil weights for the micro-vertex on that macro-edge,
///         weights are sorted according to the vertexdof-macro-edge stencil index function
///
template < class P1Form >
inline std::vector< real_t > assembleP1LocalStencil_new( const std::shared_ptr< PrimitiveStorage >& storage,
                                                         const Edge&                                edge,
                                                         const indexing::Index&                     microVertexIndex,
                                                         const uint_t&                              level,
                                                         P1Form&                                    form )
{
   // check if index lies in the edges's interior
   WALBERLA_CHECK_EQUAL( microVertexIndex.y(), 0, "[P1 edge stencil assembly] y-coordinate on edge must be zero" );
   WALBERLA_CHECK_EQUAL( microVertexIndex.z(), 0, "[P1 edge stencil assembly] z-coordinate on edge must be zero" );
   WALBERLA_CHECK_GREATER( microVertexIndex.x(), 0 );
   WALBERLA_CHECK_LESS( microVertexIndex.x(), levelinfo::num_microvertices_per_edge( level ) );

   const uint_t          stencilSize = vertexDoFMacroEdgeStencilMemorySize( level, edge );
   std::vector< real_t > stencil( stencilSize, real_c( 0 ) );

   if ( level == 0 )
      return stencil;

   for ( const auto& macroCellID : edge.neighborCells() )
   {
      const auto macroCell = storage->getCell( macroCellID );

      // 1. translate coordinate to macro-cell

      // find out the local ID of the edge in the cell
      const uint_t localEdgeID = macroCell->getLocalEdgeID( edge.getID() );

      // Find out the coordinate system basis of the index on the macro-cell.
      WALBERLA_ASSERT_EQUAL( macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].size(), 2 );
      const uint_t basisCenter     = macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 0 );
      const uint_t basisXDirection = macroCell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 1 );
      // find out the missing Z direction
      const std::set< uint_t > allDirections      = { 0, 1, 2, 3 };
      const std::set< uint_t > allDirectionsButYZ = { basisCenter, basisXDirection };
      std::vector< uint_t >    missingDirections;
      std::set_difference( allDirections.begin(),
                           allDirections.end(),
                           allDirectionsButYZ.begin(),
                           allDirectionsButYZ.end(),
                           std::inserter( missingDirections, missingDirections.begin() ) );
      WALBERLA_ASSERT_EQUAL( missingDirections.size(), 2 );
      const uint_t                  basisYDirection  = missingDirections[0];
      const uint_t                  basisZDirection  = missingDirections[1];
      const std::array< uint_t, 4 > indexingBasis    = { basisCenter, basisXDirection, basisYDirection, basisZDirection };
      const auto                    indexInMacroCell = indexing::basisConversion(
          microVertexIndex, indexingBasis, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );

      WALBERLA_DEBUG_SECTION()
      {
         const auto debugLocalEdges = vertexdof::macrocell::isOnCellEdge( indexInMacroCell, level );
         const auto debugLocalFaces = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         WALBERLA_ASSERT_EQUAL( debugLocalEdges.size(), 1 );
         WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 2 );
         WALBERLA_ASSERT_EQUAL( *debugLocalEdges.begin(), localEdgeID );
      }

      // 2. calculate stiffness matrix for each micro-cell and store contributions
      const auto cellLocalStencilWeights = calculateStencilInMacroCellForm_new( indexInMacroCell, *macroCell, level, form );

      // 3. translate coordinates / stencil directions back to edge-local coordinate system
      for ( const auto it : cellLocalStencilWeights )
      {
         const auto cellLocalDir  = it.first;
         const auto stencilWeight = it.second;

         const hyteg::indexing::Index cellLocalIndexInDir =
             indexInMacroCell + vertexdof::logicalIndexOffsetFromVertex( cellLocalDir );
         const auto            onLocalFacesCenter = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         const auto            onLocalFacesDir    = vertexdof::macrocell::isOnCellFace( cellLocalIndexInDir, level );
         std::vector< uint_t > intersectingFaces;
         std::set_intersection( onLocalFacesCenter.begin(),
                                onLocalFacesCenter.end(),
                                onLocalFacesDir.begin(),
                                onLocalFacesDir.end(),
                                std::back_inserter( intersectingFaces ) );

         if ( intersectingFaces.size() >= 2 )
         {
            // on edge
            const auto edgeLocalIndexInDir = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, indexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            WALBERLA_ASSERT_EQUAL( edgeLocalIndexInDir.y(), 0 );
            WALBERLA_ASSERT_EQUAL( edgeLocalIndexInDir.z(), 0 );
            const int dirDIfference = static_cast< int >( edgeLocalIndexInDir.x() - microVertexIndex.x() );
            WALBERLA_ASSERT_GREATER_EQUAL( dirDIfference, -1 );
            WALBERLA_ASSERT_LESS_EQUAL( dirDIfference, 1 );
            const stencilDirection dirOnEdge =
                dirDIfference == 0 ? sd::VERTEX_C : ( dirDIfference == 1 ? sd::VERTEX_E : sd::VERTEX_W );
            stencil[vertexdof::macroedge::stencilIndexOnEdge( dirOnEdge )] += stencilWeight;
         }
         else if ( intersectingFaces.size() == 1 )
         {
            // on neighbor face
            const auto localFaceIDInCell = *intersectingFaces.begin();
            const auto facePrimitiveID   = macroCell->neighborFaces()[localFaceIDInCell];
            // To get the correct indexing basis, we check which one results in a zero entry in the z coordinate.
            const auto                    firstTestIndexingBasis  = indexingBasis;
            const std::array< uint_t, 4 > secondTestIndexingBasis = {
                indexingBasis[0], indexingBasis[1], indexingBasis[3], indexingBasis[2] };
            const auto edgeLocalIndexInDirFirst = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, firstTestIndexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            const auto edgeLocalIndexInDirSecond = indexing::basisConversion(
                cellLocalIndexInDir, { 0, 1, 2, 3 }, secondTestIndexingBasis, levelinfo::num_microvertices_per_edge( level ) );
            WALBERLA_ASSERT_UNEQUAL( edgeLocalIndexInDirFirst.z(), edgeLocalIndexInDirSecond.z() );
            WALBERLA_ASSERT( edgeLocalIndexInDirFirst.z() == 0 || edgeLocalIndexInDirSecond.z() == 0 );
            const hyteg::indexing::Index faceLocalIndexInDir =
                edgeLocalIndexInDirFirst.z() == 0 ? edgeLocalIndexInDirFirst : edgeLocalIndexInDirSecond;
            WALBERLA_ASSERT_EQUAL( faceLocalIndexInDir.y(), 1 );
            stencilDirection faceLocalStencilDirection;

            const auto xOffset = static_cast< int >( faceLocalIndexInDir.x() ) - static_cast< int >( microVertexIndex.x() );
            if ( xOffset == 0 )
            {
               faceLocalStencilDirection = sd::VERTEX_E;
            }
            else if ( xOffset == -1 )
            {
               faceLocalStencilDirection = sd::VERTEX_W;
            }
            else
            {
               WALBERLA_ABORT( "[P1Elements][Edge] Invalid offsets" );
            }
            stencil[vertexdof::macroedge::stencilIndexOnNeighborFace( faceLocalStencilDirection,
                                                                      edge.face_index( facePrimitiveID ) )] += stencilWeight;
         }
         else if ( intersectingFaces.size() == 0 )
         {
            // in macro-cell
            stencil[vertexdof::macroedge::stencilIndexOnNeighborCell( edge.cell_index( macroCellID ),
                                                                      edge.getNumNeighborFaces() )] += stencilWeight;
         }
      }
   }
   return stencil;
}

/// \brief Assembles the local P1 operator stencil on a macro-face
///
/// \param storage the governing \ref PrimitiveStorage
/// \param face the macro-face
/// \param microVertexIndex the micro-vertex index on the macro-face (z coordinate must be 0)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a map containing the stencil weights for the micro-vertex on that macro-face,
///         it maps the 13 or 19 (== 7 on face + 6 on ghost face 0 + 6 on ghost face 1) stencil directions to the stencil weights
///         refer to the documentation to better understand that mapping
///
template < class P1Form >
inline std::map< stencilDirection, real_t > assembleP1LocalStencil( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    const Face&                                face,
                                                                    const indexing::Index&                     microVertexIndex,
                                                                    const uint_t&                              level,
                                                                    const P1Form&                              form )
{
   // check if index lies in the face's interior
   WALBERLA_CHECK_EQUAL( microVertexIndex.z(), 0, "[P1 face stencil assembly] z-coordinate on face must be zero" );
   WALBERLA_CHECK_GREATER( microVertexIndex.x() + microVertexIndex.y(), 0 );
   WALBERLA_CHECK_LESS( microVertexIndex.x() + microVertexIndex.y(), levelinfo::num_microvertices_per_edge( level ) );

   std::map< stencilDirection, real_t > faceStencil;

   for ( const auto& macroCellID : face.neighborCells() )
   {
      const auto macroCell = storage->getCell( macroCellID );

      // 1. translate coordinate to macro-cell

      // find out the local ID of the face in the cell
      const uint_t localFaceID = macroCell->getLocalFaceID( face.getID() );

      // find out the coordinate system basis of the index on the macro-cell
      WALBERLA_ASSERT_EQUAL( macroCell->getFaceLocalVertexToCellLocalVertexMaps()[localFaceID].size(), 3 );
      const uint_t basisCenter     = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[localFaceID].at( 0 );
      const uint_t basisXDirection = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[localFaceID].at( 1 );
      const uint_t basisYDirection = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[localFaceID].at( 2 );
      // find out the missing Z direction
      const std::set< uint_t > allDirections     = { 0, 1, 2, 3 };
      const std::set< uint_t > allDirectionsButZ = { basisCenter, basisXDirection, basisYDirection };
      std::set< uint_t >       missingDirection;
      std::set_difference( allDirections.begin(),
                           allDirections.end(),
                           allDirectionsButZ.begin(),
                           allDirectionsButZ.end(),
                           std::inserter( missingDirection, missingDirection.begin() ) );
      WALBERLA_ASSERT_EQUAL( missingDirection.size(), 1 );
      const uint_t                  basisZDirection  = *missingDirection.begin();
      const std::array< uint_t, 4 > indexingBasis    = { basisCenter, basisXDirection, basisYDirection, basisZDirection };
      const auto                    indexInMacroCell = indexing::basisConversion(
          microVertexIndex, indexingBasis, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );

      WALBERLA_DEBUG_SECTION()
      {
         const auto debugLocalFaces = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
         WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 1 );
         WALBERLA_ASSERT_EQUAL( *debugLocalFaces.begin(), localFaceID );
      }

      // 2. calculate stiffness matrix for each micro-cell and store contributions
      const auto cellLocalStencilWeights = calculateStencilInMacroCellForm( indexInMacroCell, *macroCell, level, form );

      // 3. translate coordinates / stencil directions back to face-local coordinate system
      for ( const auto it : cellLocalStencilWeights )
      {
         const auto cellLocalDir  = it.first;
         const auto stencilWeight = it.second;

         const auto cellLocalIndexInDir = indexInMacroCell + vertexdof::logicalIndexOffsetFromVertex( cellLocalDir );
         const auto faceLocalIndexInDir = indexing::basisConversion(
             cellLocalIndexInDir, { 0, 1, 2, 3 }, indexingBasis, levelinfo::num_microvertices_per_edge( level ) );
         WALBERLA_ASSERT_LESS_EQUAL( faceLocalIndexInDir.z(), 1 );
         const auto indexOnGhostLayer = faceLocalIndexInDir.z() == 1;
         const auto localCellID       = face.cell_index( macroCellID );
         WALBERLA_ASSERT_LESS_EQUAL( localCellID, 1 );
         const auto faceLocalStencilDirection = [&face, microVertexIndex, faceLocalIndexInDir, indexOnGhostLayer, localCellID] {
            const auto       xOffset = static_cast< int >( faceLocalIndexInDir.x() ) - static_cast< int >( microVertexIndex.x() );
            const auto       yOffset = static_cast< int >( faceLocalIndexInDir.y() ) - static_cast< int >( microVertexIndex.y() );
            stencilDirection projectedDirection;
            if ( xOffset == 0 && yOffset == 0 )
               projectedDirection = stencilDirection::VERTEX_C;
            else if ( xOffset == 1 && yOffset == 1 )
               projectedDirection = stencilDirection::VERTEX_NE;
            else if ( xOffset == 0 && yOffset == 1 )
               projectedDirection = stencilDirection::VERTEX_N;
            else if ( xOffset == -1 && yOffset == 1 )
               projectedDirection = stencilDirection::VERTEX_NW;
            else if ( xOffset == 1 && yOffset == 0 )
               projectedDirection = stencilDirection::VERTEX_E;
            else if ( xOffset == -1 && yOffset == 0 )
               projectedDirection = stencilDirection::VERTEX_W;
            else if ( xOffset == 1 && yOffset == -1 )
               projectedDirection = stencilDirection::VERTEX_SE;
            else if ( xOffset == 0 && yOffset == -1 )
               projectedDirection = stencilDirection::VERTEX_S;
            else if ( xOffset == -1 && yOffset == -1 )
               projectedDirection = stencilDirection::VERTEX_SW;
            else
            {
               WALBERLA_ASSERT( false, "Invalid offsets" );
               projectedDirection = stencilDirection::VERTEX_TC;
            }

            if ( indexOnGhostLayer )
            {
               // deciding here that the stencil direction for the first cell at a face is top
               if ( localCellID == 0 )
                  return makeVertexDirectionTop( projectedDirection );
               else
               {
                  WALBERLA_ASSERT_EQUAL( face.getNumNeighborCells(), 2 );
                  WALBERLA_UNUSED( face );
                  return makeVertexDirectionBottom( projectedDirection );
               }
            }
            else
            {
               return projectedDirection;
            }
         }();

         if ( faceStencil.count( faceLocalStencilDirection ) == 0 )
         {
            faceStencil[faceLocalStencilDirection] = real_c( 0 );
         }
         faceStencil[faceLocalStencilDirection] += stencilWeight;
      }
   }
   return faceStencil;
}

/// \brief Assembles the local P1 operator stencil on a macro-cell
///
/// \param storage the governing \ref PrimitiveStorage
/// \param cell the macro-cell
/// \param microVertexIndex the micro-vertex index on the macro-cell (must lie in the interior of the macro-cell)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a map containing the stencil weights for the micro-vertex on that macro-cell,
///
template < class P1Form >
inline std::map< stencilDirection, real_t > assembleP1LocalStencil( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    const Cell&                                cell,
                                                                    const indexing::Index&                     microVertexIndex,
                                                                    const uint_t&                              level,
                                                                    const P1Form&                              form )
{
   WALBERLA_UNUSED( storage );
   WALBERLA_DEBUG_SECTION()
   {
      const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
      const auto onCellEdges    = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
      const auto onCellFaces    = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );
      WALBERLA_CHECK_EQUAL( onCellVertices.size(), 0 );
      WALBERLA_CHECK_EQUAL( onCellEdges.size(), 0 );
      WALBERLA_CHECK_EQUAL( onCellFaces.size(), 0 );
   }
   return calculateStencilInMacroCellForm( microVertexIndex, cell, level, form );
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
/// \brief Assembles the local P1 operator stencil on a macro-cell
///
/// \param storage the governing \ref PrimitiveStorage
/// \param cell the macro-cell
/// \param microVertexIndex the micro-vertex index on the macro-cell (must lie in the interior of the macro-cell)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a map containing the stencil weights for the micro-vertex on that macro-cell,
///
template < class P1Form >
inline std::map< stencilDirection, real_t > assembleP1LocalStencil_new( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        const Cell&                                cell,
                                                                        const indexing::Index& microVertexIndex,
                                                                        const uint_t&          level,
                                                                        P1Form&                form )
{
   WALBERLA_UNUSED( storage );
   WALBERLA_DEBUG_SECTION()
   {
      const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
      const auto onCellEdges    = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
      const auto onCellFaces    = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );
      WALBERLA_CHECK_EQUAL( onCellVertices.size(), 0 );
      WALBERLA_CHECK_EQUAL( onCellEdges.size(), 0 );
      WALBERLA_CHECK_EQUAL( onCellFaces.size(), 0 );
   }
   return calculateStencilInMacroCellForm_new( microVertexIndex, cell, level, form );
}

/// \brief Assembles the local P1 operator stencil on a macro-cell
///
/// \param storage the governing \ref PrimitiveStorage
/// \param cell the macro-cell
/// \param microVertexIndex the micro-vertex index on the macro-cell (must lie in the interior of the macro-cell)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a map containing the stencil weights for the micro-vertex on that macro-cell,
///
template < class P1Form >
inline std::map< indexing::Index, real_t > assembleP1LocalStencilNew( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                               const Cell&                                cell,
                                                                               const indexing::Index& microVertexIndex,
                                                                               const uint_t&          level,
                                                                               const P1Form&          form )
{
   WALBERLA_UNUSED( storage );
   WALBERLA_DEBUG_SECTION()
   {
      //    const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
      //    const auto onCellEdges = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
      //    const auto onCellFaces = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );
      //    WALBERLA_CHECK_EQUAL( onCellVertices.size(), 0 );
      //    WALBERLA_CHECK_EQUAL( onCellEdges.size(), 0 );
      //    WALBERLA_CHECK_EQUAL( onCellFaces.size(), 0 );
   }
   auto stencilMap = calculateStencilInMacroCellForm( microVertexIndex, cell, level, form );
   std::map< indexing::Index, real_t > convertedMap;
   for ( const auto& it : stencilMap )
   {
      convertedMap[vertexdof::logicalIndexOffsetFromVertex( it.first )] = it.second;
   }
   return convertedMap;
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
/// \brief Assembles the local P1 operator stencil on a macro-cell
///
/// \param storage the governing \ref PrimitiveStorage
/// \param cell the macro-cell
/// \param microVertexIndex the micro-vertex index on the macro-cell (must lie in the interior of the macro-cell)
/// \param level the multigrid level
/// \param form the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
/// \return a map containing the stencil weights for the micro-vertex on that macro-cell,
///
template < class P1Form >
inline std::map< indexing::Index, real_t >
    assembleP1LocalStencilNew_new( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const Cell&                                cell,
                                   const indexing::Index&                     microVertexIndex,
                                   const uint_t&                              level,
                                   P1Form&                                    form )
{
   WALBERLA_UNUSED( storage );
   WALBERLA_DEBUG_SECTION()
   {
      //    const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
      //    const auto onCellEdges = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
      //    const auto onCellFaces = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );
      //    WALBERLA_CHECK_EQUAL( onCellVertices.size(), 0 );
      //    WALBERLA_CHECK_EQUAL( onCellEdges.size(), 0 );
      //    WALBERLA_CHECK_EQUAL( onCellFaces.size(), 0 );
   }
   auto stencilMap = calculateStencilInMacroCellForm_new( microVertexIndex, cell, level, form );
   std::map< indexing::Index, real_t > convertedMap;
   for ( const auto& it : stencilMap )
   {
      convertedMap[vertexdof::logicalIndexOffsetFromVertex( it.first )] = it.second;
   }
   return convertedMap;
}

} // namespace P1Elements3D
} // namespace P1Elements
} // namespace hyteg

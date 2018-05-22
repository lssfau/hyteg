#pragma once

#include <unordered_map>
#include <tinyhhg_core/celldofspace/CellDoFIndexing.hpp>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"


namespace hhg {

using walberla::real_t;
using walberla::real_c;
using walberla::uint_t;

namespace P1Elements {
namespace FaceVertexDoF {

// Fenics P1 DoF ordering
// 2       1---0
// |\       \  |
// | \       \ |
// |  \       \|
// 0---1       2

const uint_t ElementSize = 3;

typedef stencilDirection SD;
typedef std::array<SD, ElementSize> P1Element;
typedef std::array<uint_t, ElementSize> DoFMap;
typedef std::array<uint_t, ElementSize> StencilMap;

const P1Element elementSW = {{SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S}};
const P1Element elementS = {{SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE}};
const P1Element elementSE = {{SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E}};
const P1Element elementNE = {{SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N}};
const P1Element elementN = {{SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_NW}};
const P1Element elementNW = {{SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_W}};

// ordered
const P1Element elementSWOrd = {{SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S}};
const P1Element elementSOrd = {{SD::VERTEX_S, SD::VERTEX_SE, SD::VERTEX_C}};
const P1Element elementSEOrd = {{SD::VERTEX_E, SD::VERTEX_C, SD::VERTEX_SE}};
const P1Element elementNEOrd = {{SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N}};
const P1Element elementNOrd = {{SD::VERTEX_N, SD::VERTEX_NW, SD::VERTEX_C}};
const P1Element elementNWOrd = {{SD::VERTEX_W, SD::VERTEX_C, SD::VERTEX_NW}};

static const std::array<P1Element, 3> P1GrayElements =
    {{
         elementS,
         elementNE,
         elementNW
     }};

static const std::array<P1Element, 3> P1BlueElements =
    {{
         elementSW,
         elementSE,
         elementN
     }};

static const std::array<StencilMap, 3> P1GrayStencilMaps =
    {{
         {{vertexdof::stencilIndexFromVertex(elementS[0]), vertexdof::stencilIndexFromVertex(elementS[1]), vertexdof::stencilIndexFromVertex(elementS[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNE[0]), vertexdof::stencilIndexFromVertex(elementNE[1]), vertexdof::stencilIndexFromVertex(elementNE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNW[0]), vertexdof::stencilIndexFromVertex(elementNW[1]), vertexdof::stencilIndexFromVertex(elementNW[2])}}
     }};

static const std::array<StencilMap, 3> P1BlueStencilMaps =
    {{
         {{vertexdof::stencilIndexFromVertex(elementSW[0]), vertexdof::stencilIndexFromVertex(elementSW[1]), vertexdof::stencilIndexFromVertex(elementSW[2])}},
         {{vertexdof::stencilIndexFromVertex(elementSE[0]), vertexdof::stencilIndexFromVertex(elementSE[1]), vertexdof::stencilIndexFromVertex(elementSE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementN[0]), vertexdof::stencilIndexFromVertex(elementN[1]), vertexdof::stencilIndexFromVertex(elementN[2])}}
     }};

static const std::array<DoFMap, 3> P1GrayDoFMaps =
    {{
         {{2, 0, 1}},
         {{0, 1, 2}},
         {{1, 2, 0}}
     }};


static const std::array<DoFMap, 3> P1BlueDoFMaps =
    {{
         {{0, 1, 2}},
         {{1, 2, 0}},
         {{2, 0, 1}}
     }};

inline StencilMap convertStencilDirectionsToIndices( const P1Element & element )
{
  return {{ vertexdof::stencilIndexFromVertex( element[0] ), vertexdof::stencilIndexFromVertex( element[1] ), vertexdof::stencilIndexFromVertex( element[2] ) }};
}

template<typename StencilMemory>
inline void assembleP1LocalStencil(const StencilMap &stencilMap, const DoFMap &dofMap, const Matrix3r &localMatrix,
                                   StencilMemory &stencil, double coeffWeight = 1.0) {
  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += coeffWeight * localMatrix(dofMap[0], dofMap[j]);
  }
}

}





namespace CellVertexDoF {

typedef stencilDirection sd;

const std::array< std::array< stencilDirection, 4 >, 4 > whiteUpCellsAtInnerVertex = {{
                                                                                 { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN }, // below
                                                                                 { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS }, // top front
                                                                                 { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW }, // top back west
                                                                                 { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC }, // top back east
                                                                                 }};

const std::array< std::array< stencilDirection, 4 >, 4 > whiteDownCellsAtInnerVertex = {{
                                                                                   { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_S }, // below front west
                                                                                   { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_BE }, // below front east
                                                                                   { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_NW, sd::VERTEX_BN }, // below back
                                                                                   { sd::VERTEX_C, sd::VERTEX_TS, sd::VERTEX_TC, sd::VERTEX_TW }, // top
                                                                                   }};

const std::array< std::array< stencilDirection, 4 >, 4 > blueUpCellsAtInnerVertex = {{
                                                                                { sd::VERTEX_C, sd::VERTEX_BC,  sd::VERTEX_BN, sd::VERTEX_BNW  }, // below
                                                                                { sd::VERTEX_C, sd::VERTEX_W,  sd::VERTEX_S, sd::VERTEX_TS  }, // top front west
                                                                                { sd::VERTEX_C, sd::VERTEX_E,  sd::VERTEX_SE, sd::VERTEX_TSE }, // top front east
                                                                                { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC   }, // top back
                                                                                }};

const std::array< std::array< stencilDirection, 4 >, 4 > blueDownCellsAtInnerVertex = {{
                                                                                  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE }, // below front
                                                                                  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW }, // below back west
                                                                                  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N }, // below back east
                                                                                  { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE }, // top
                                                                                  }};

const std::array< std::array< stencilDirection, 4 >, 4 > greenUpCellsAtInnerVertex = {{
                                                                                 { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW }, // below west
                                                                                 { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN }, // below east
                                                                                 { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW }, // top back
                                                                                 { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE }, // top front
                                                                                 }};

const std::array< std::array< stencilDirection, 4 >, 4 > greenDownCellsAtInnerVertex = {{
                                                                                   { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE }, // below front
                                                                                   { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW }, // below back
                                                                                   { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC }, // top east
                                                                                   { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW }, // top west
                                                                                   }};


const std::array< std::array< stencilDirection, 4 >, 24 > allCellsAtInnerVertex = {{
                                                                        whiteUpCellsAtInnerVertex[0], whiteUpCellsAtInnerVertex[1], whiteUpCellsAtInnerVertex[2], whiteUpCellsAtInnerVertex[3],
                                                                        whiteDownCellsAtInnerVertex[0], whiteDownCellsAtInnerVertex[1], whiteDownCellsAtInnerVertex[2], whiteDownCellsAtInnerVertex[3],
                                                                        blueUpCellsAtInnerVertex[0], blueUpCellsAtInnerVertex[1], blueUpCellsAtInnerVertex[2], blueUpCellsAtInnerVertex[3],
                                                                        blueDownCellsAtInnerVertex[0], blueDownCellsAtInnerVertex[1], blueDownCellsAtInnerVertex[2], blueDownCellsAtInnerVertex[3],
                                                                        greenUpCellsAtInnerVertex[0], greenUpCellsAtInnerVertex[1], greenUpCellsAtInnerVertex[2], greenUpCellsAtInnerVertex[3],
                                                                        greenDownCellsAtInnerVertex[0], greenDownCellsAtInnerVertex[1], greenDownCellsAtInnerVertex[2], greenDownCellsAtInnerVertex[3]
                                                                        }};


const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace0 = {{
  // no cells with bottom direction
  { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
  { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
  { sd::VERTEX_C, sd::VERTEX_TS, sd::VERTEX_TC, sd::VERTEX_TW },
  { sd::VERTEX_C, sd::VERTEX_W,  sd::VERTEX_S, sd::VERTEX_TS  },
  { sd::VERTEX_C, sd::VERTEX_E,  sd::VERTEX_SE, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC   },
  { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW },
  { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW }
}};

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace1 = {{
  // no cells with south direction
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
  { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
  { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_NW, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_BC,  sd::VERTEX_BN, sd::VERTEX_BNW  },
  { sd::VERTEX_C, sd::VERTEX_NW, sd::VERTEX_N, sd::VERTEX_TC   },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TW, sd::VERTEX_NW },
  { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW }
}};

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace2 = {{
  // no cells with west direction
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
  { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_E, sd::VERTEX_TC },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_BE },
  { sd::VERTEX_C, sd::VERTEX_E,  sd::VERTEX_SE, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BN, sd::VERTEX_N },
  { sd::VERTEX_C, sd::VERTEX_TC, sd::VERTEX_TS, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BE, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_SE, sd::VERTEX_TS, sd::VERTEX_TSE },
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE },
  { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_TSE, sd::VERTEX_TC }
}};

const std::array< std::array< stencilDirection, 4 >, 12 > allCellsAtFace3 = {{
  // no cells in {N, E, TC, TSE}
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_BN },
  { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_TS },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_TW },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_S },
  { sd::VERTEX_C, sd::VERTEX_BC,  sd::VERTEX_BN, sd::VERTEX_BNW  },
  { sd::VERTEX_C, sd::VERTEX_W,  sd::VERTEX_S, sd::VERTEX_TS  },
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_S, sd::VERTEX_SE },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_NW, sd::VERTEX_BNW },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BC, sd::VERTEX_BNW },
  { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BE, sd::VERTEX_SE },
  { sd::VERTEX_C, sd::VERTEX_BN, sd::VERTEX_BNW, sd::VERTEX_NW },
  { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_TS, sd::VERTEX_TW },
}};


inline std::vector< std::array< stencilDirection, 4 > > getNeighboringElements( const indexing::Index & microVertexIndex, const uint_t & level )
{
  typedef std::vector< std::array< stencilDirection, 4 > > returnType;

  const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( microVertexIndex, level );
  const auto onCellEdges    = vertexdof::macrocell::isOnCellEdge( microVertexIndex, level );
  const auto onCellFaces    = vertexdof::macrocell::isOnCellFace( microVertexIndex, level );

  if ( onCellVertices.size() > 0 )
  {
    WALBERLA_ABORT( "Stencil assembly (3D) on macro-vertices not implemented" );
  }
  else if ( onCellEdges.size() > 0 )
  {
    WALBERLA_ABORT( "Stencil assembly (3D) on macro-edges not implemented" );
  }
  else if ( onCellFaces.size() > 0 )
  {
    WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 1 );
    const auto localFaceID = *onCellFaces.begin();
    WALBERLA_ASSERT_GREATER_EQUAL( localFaceID, 0 );
    WALBERLA_ASSERT_LESS_EQUAL   ( localFaceID, 3 );
    switch( localFaceID )
    {
      case 0:
        return returnType( allCellsAtFace0.begin(), allCellsAtFace0.end() );
      case 1:
        return returnType( allCellsAtFace1.begin(), allCellsAtFace1.end() );
      case 2:
        return returnType( allCellsAtFace2.begin(), allCellsAtFace2.end() );
      case 3:
        return returnType( allCellsAtFace3.begin(), allCellsAtFace3.end() );
      default:
        return returnType( allCellsAtFace0.begin(), allCellsAtFace0.end() );
    }
  }
  else
  {
    return returnType( allCellsAtInnerVertex.begin(), allCellsAtInnerVertex.end() );
  }
};


/// Calculates the stencil weights from the stiffness matrices of neighboring elements at an index in a macro-cell.
/// Also works for indices on the boundary of a macro-cell.
/// It automatically computes / selects the neighboring elements depending on the micro-vertex' location.
/// \param microVertexIndex the logical index of the micro-vertex in a macro-cell (can alos be located on the macro-cell's boundary)
/// \param cell the surrounding macro-cell
/// \param level the hierarchy level
/// \param ufcGen the UFC object that implements tabulate_tensor() to calculate the local stiffness matrix
template< typename UFCOperator >
inline std::map< stencilDirection, real_t > calculateStencilInMacroCell( const indexing::Index & microVertexIndex, const Cell & cell,
                                                                         const uint_t & level, const UFCOperator & ufcGen )
{
  std::map< stencilDirection, real_t > macroCellStencilEntries;

  const auto neighboringElements = getNeighboringElements( microVertexIndex, level );

  // 1. Going over all neighboring cells of a micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
  //    VERTEX_C represents the reference micro-vertex.
  for ( const auto & cellAtVertex : neighboringElements )
  {
    WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

    // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
    std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ ) {
      logicalOffsetsFromCenter[localID] = microVertexIndex + vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[localID] );
    }

    // 3. Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
    std::array< Point3D, 4 > geometricCoordinates;
    for ( uint_t localID = 0; localID < 4; localID++ ) {
      geometricCoordinates[localID] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[localID] );
    }

    std::array< Point3D, 4 > geometricOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ ) {
      geometricOffsetsFromCenter[localID] = geometricCoordinates[localID] - geometricCoordinates[0];
    }

    // 4. Computing the local stiffness matrix
    //    To calculate the 4x4 stiffness matrix, we need the geometric offsets from the reference micro-vertex
    //    from all micro-vertices in the neighbor cell (including the reference micro-vertex itself -> the first offset is always (0.0, 0.0, 0.0))

    // Flattening the offset array to be able to pass it to the fenics routines.
    double geometricOffsetsArray[12];
    for ( uint_t cellVertex = 0; cellVertex < 4; cellVertex++ ) {
      for ( uint_t coordinate = 0; coordinate < 3; coordinate++ ) {
        geometricOffsetsArray[cellVertex * 3 + coordinate] = geometricOffsetsFromCenter[cellVertex][coordinate];
      }
    }

    Matrix4r localStiffnessMatrix;
    ufcGen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

    // 5. Adding contribution to stencil
    //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
    //    we always only need the first row of the local stiffness matrix.
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      const stencilDirection stencilDir = cellAtVertex[ localID ];
      if ( macroCellStencilEntries.count( stencilDir ) == 0 )
      {
        macroCellStencilEntries[ stencilDir ] = real_c( 0 );
      }
      macroCellStencilEntries[ stencilDir ] += localStiffnessMatrix( 0, localID );
    }
  }
  return macroCellStencilEntries;
}


template< typename UFCOperator >
inline std::array< real_t, 15 /* is 15 enough? */ > assembleP1LocalStencil( const std::shared_ptr< PrimitiveStorage > & storage, const Face & face,
                                                                            const indexing::Index & microVertexIndex, const uint_t & level, const UFCOperator & ufcGen )
{
  // TODO: check if index in the face's interior

  for ( const auto & macroCellID : face.neighborCells() )
  {
    const auto macroCell = storage->getCell( macroCellID );

    // 1. translate coordinate to macro-cell
    
    // find out the local ID of the face in the cell
    const uint_t localFaceID = macroCell->getLocalFaceID( face.getID() );
    
    // find out the coordinate system basis of the index on the macro-cell
    WALBERLA_ASSERT_EQUAL( macroCell->getFaceLocalVertexToCellLocalVertexMaps()[ localFaceID ].size(), 3 );
    const uint_t basisCenter     = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[ localFaceID ].at( 0 );
    const uint_t basisXDirection = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[ localFaceID ].at( 1 );
    const uint_t basisYDirection = macroCell->getFaceLocalVertexToCellLocalVertexMaps()[ localFaceID ].at( 2 );
    // find out the missing Z direction
    const std::set< uint_t > allDirections     = { 0, 1, 2, 3 };
    const std::set< uint_t > allDirectionsButZ = { basisCenter, basisXDirection, basisYDirection };
    std::set< uint_t > missingDirection;
    std::set_difference( allDirections.begin(), allDirections.end(), allDirectionsButZ.begin(), allDirectionsButZ.end(), std::inserter(missingDirection, missingDirection.begin()) );
    WALBERLA_ASSERT_EQUAL( missingDirection.size(), 1 );
    const uint_t basisZDirection = *missingDirection.begin();
    const std::array< uint_t, 4 > indexingBasis = { basisCenter, basisXDirection, basisYDirection, basisZDirection };
    const auto indexInMacroCell = indexing::basisConversion( microVertexIndex, { 0, 1, 2, 3 }, indexingBasis, level );

    WALBERLA_DEBUG_SECTION()
    {
      const auto debugLocalFaces = vertexdof::macrocell::isOnCellFace( indexInMacroCell, level );
      WALBERLA_ASSERT_EQUAL( debugLocalFaces.size(), 1 );
      WALBERLA_ASSERT_EQUAL( *debugLocalFaces.begin(), localFaceID );
    }

    // 2. find neighboring micro-cells
    const auto cellsAtFace = [ localFaceID ]
    {
        WALBERLA_ASSERT_GREATER_EQUAL( localFaceID, 0, "Invalid local face ID." );
        WALBERLA_ASSERT_LESS_EQUAL   ( localFaceID, 3, "Invalid local face ID." );
        switch ( localFaceID )
        {
          case 0:
            return allCellsAtFace0;
          case 1:
            return allCellsAtFace1;
          case 2:
            return allCellsAtFace2;
          case 3:
            return allCellsAtFace3;
          default:
            return allCellsAtFace0;
        }
    }();

    // 3. calculate stiffness matrix for each micro-cell and store contributions
    const auto cellLocalStencilWeights = calculateStencilInMacroCell( indexInMacroCell, *macroCell, level, ufcGen );

    // 4. TODO: translate coordinates / stencil directions back to face-local coordinate system

  }

}

/// Assembles the (constant) P1 stencil for one macro-cell using the local stiffness matrix calculated by the passed fenics UFC generator.
/// \param cell the macro-cell for which the stencil shall be assembled
/// \param level the HHG refinement level
/// \param ufcGen a cell_integral subclass from the generated output of the fenics library
///               which implements the member tabulate_tensor() that calculates the local stiffness matrix
/// \return an array of stencil weights - the stencil weights are sorted according to the stencil index functions
template< typename UFCOperator >
inline std::array< real_t, 15 > assembleP1LocalStencil( const Cell & cell, const uint_t & level, const UFCOperator & ufcGen )
{
  std::array< real_t, 15 > localStencil;
  localStencil.fill( real_c( 0 ) );

  const auto stencilMap = calculateStencilInMacroCell( indexing::Index( 1, 1, 1 ), cell, level, ufcGen );
  for ( const auto stencilEntry : stencilMap )
  {
    const auto dir    = stencilEntry.first;
    const auto weight = stencilEntry.second;
    const uint_t stencilIndex = vertexdof::stencilIndexFromVertex( dir );
    localStencil[ stencilIndex ] = weight;
  }

  return localStencil;
}

}
}
}

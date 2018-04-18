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

namespace P1Elements {

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

namespace FaceVertexDoF {

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
}

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

namespace CellVertexDoF {

typedef stencilDirection sd;

const std::array< std::array< stencilDirection, 4 >, 4 > whiteUpCellsAtVertex = {{
                                                                           { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_BS }, // below
                                                                           { sd::VERTEX_C, sd::VERTEX_FC, sd::VERTEX_FE, sd::VERTEX_FN }, // top front
                                                                           { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BW, sd::VERTEX_NW }, // top back west
                                                                           { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_E, sd::VERTEX_N }, // top back east
                                                                           }};

const std::array< std::array< stencilDirection, 4 >, 4 > whiteDownCellsAtVertex = {{
                                                                             { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_S, sd::VERTEX_FC }, // below front west
                                                                             { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_FE, sd::VERTEX_SE }, // below front east
                                                                             { sd::VERTEX_C, sd::VERTEX_BC, sd::VERTEX_BW, sd::VERTEX_BS }, // below back
                                                                             { sd::VERTEX_C, sd::VERTEX_FN, sd::VERTEX_N, sd::VERTEX_NW }, // top
                                                                             }};

const std::array< std::array< stencilDirection, 4 >, 4 > blueUpCellsAtVertex = {{
                                                                          { sd::VERTEX_C, sd::VERTEX_S,  sd::VERTEX_BS, sd::VERTEX_BSW  }, // below
                                                                          { sd::VERTEX_C, sd::VERTEX_W,  sd::VERTEX_FC, sd::VERTEX_FN  }, // top front west
                                                                          { sd::VERTEX_C, sd::VERTEX_E,  sd::VERTEX_FE, sd::VERTEX_FNE }, // top front east
                                                                          { sd::VERTEX_C, sd::VERTEX_BW, sd::VERTEX_BC, sd::VERTEX_N   }, // top back
                                                                          }};

const std::array< std::array< stencilDirection, 4 >, 4 > blueDownCellsAtVertex = {{
                                                                            { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_FC, sd::VERTEX_FE }, // below front
                                                                            { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_BW, sd::VERTEX_BSW }, // below back west
                                                                            { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_BS, sd::VERTEX_BC }, // below back east
                                                                            { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_FN, sd::VERTEX_FNE }, // top
                                                                            }};

const std::array< std::array< stencilDirection, 4 >, 4 > greenUpCellsAtVertex = {{
                                                                           { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_S, sd::VERTEX_BSW }, // below west
                                                                           { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_SE, sd::VERTEX_BS }, // below east
                                                                           { sd::VERTEX_C, sd::VERTEX_N, sd::VERTEX_NW, sd::VERTEX_BW }, // top back
                                                                           { sd::VERTEX_C, sd::VERTEX_FE, sd::VERTEX_FN, sd::VERTEX_FNE }, // top front
                                                                           }};

const std::array< std::array< stencilDirection, 4 >, 4 > greenDownCellsAtVertex = {{
                                                                             { sd::VERTEX_C, sd::VERTEX_S, sd::VERTEX_SE, sd::VERTEX_FE }, // below front
                                                                             { sd::VERTEX_C, sd::VERTEX_BS, sd::VERTEX_BSW, sd::VERTEX_BW }, // below back
                                                                             { sd::VERTEX_C, sd::VERTEX_E, sd::VERTEX_FNE, sd::VERTEX_N }, // top east
                                                                             { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_FN, sd::VERTEX_NW }, // top west
                                                                             }};

const std::array< std::array< stencilDirection, 4 >, 24 > allCellsAtVertex = {{
                                                                        whiteUpCellsAtVertex[0], whiteUpCellsAtVertex[1], whiteUpCellsAtVertex[2], whiteUpCellsAtVertex[3],
                                                                        whiteDownCellsAtVertex[0], whiteDownCellsAtVertex[1], whiteDownCellsAtVertex[2], whiteDownCellsAtVertex[3],
                                                                        blueUpCellsAtVertex[0], blueUpCellsAtVertex[1], blueUpCellsAtVertex[2], blueUpCellsAtVertex[3],
                                                                        blueDownCellsAtVertex[0], blueDownCellsAtVertex[1], blueDownCellsAtVertex[2], blueDownCellsAtVertex[3],
                                                                        greenUpCellsAtVertex[0], greenUpCellsAtVertex[1], greenUpCellsAtVertex[2], greenUpCellsAtVertex[3],
                                                                        greenDownCellsAtVertex[0], greenDownCellsAtVertex[1], greenDownCellsAtVertex[2], greenDownCellsAtVertex[3]
                                                                        }};


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

  // 1. Going over all neighboring cells of a reference micro-vertex
  //    A neighboring cell is defined by a 4-tuple of (different) stencil directions with one of them being VERTEX_C.
  //    VERTEX_C represents the reference micro-vertex.
  for ( const auto & cellAtVertex : allCellsAtVertex )
  {
    // 2. Collecting the logical index offsets of each micro-vertex of the current neighboring cell from the reference micro-vertex
    //    The reference micro-vertex is chosen so that its logical coordinates are ( 1, 1, 1 ).
    //    This is because we do not need to use negative coordinates during the calculations and it works for all levels >= 2.
    std::array< indexing::Index, 4 > logicalOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      logicalOffsetsFromCenter[ localID ] = indexing::Index( 1, 1, 1 ) + vertexdof::logicalIndexOffsetFromVertex( cellAtVertex[ localID ] );
    }

    // 3. Calculating the absolute offsets of each micro-vertex of the current cell from the reference micro-vertex
    std::array< Point3D, 4 > geometricCoordinates;
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      geometricCoordinates[ localID ] = vertexdof::macrocell::coordinateFromIndex( level, cell, logicalOffsetsFromCenter[ localID ] );
    }

    std::array< Point3D, 4 > geometricOffsetsFromCenter;
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      WALBERLA_ASSERT_EQUAL( cellAtVertex[ 0 ], sd::VERTEX_C );
      geometricOffsetsFromCenter[ localID ] = geometricCoordinates[ localID ] - geometricCoordinates[ 0 ];
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
        geometricOffsetsArray[ cellVertex * 3 + coordinate ] = geometricOffsetsFromCenter[ cellVertex ][ coordinate ];
      }
    }

    Matrix4r localStiffnessMatrix;
    ufcGen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

    // 5. Adding contribution to stencil
    //    Since we enforced that the first entry in the local cell micro-vertex array is always the reference micro-vertex
    //    we always only need the first row of the local stiffness matrix.
    for ( uint_t localID = 0; localID < 4; localID++ )
    {
      const sd     stencilDir   = cellAtVertex[ localID ];
      const uint_t stencilIndex = vertexdof::stencilIndexFromVertex( stencilDir );

      WALBERLA_ASSERT_EQUAL( cellAtVertex[0], sd::VERTEX_C );

      localStencil[ stencilIndex ] += localStiffnessMatrix( 0, localID );
    }
  }

  return localStencil;
}

}
}
}

#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/dgfunctionspace/DGEdgeIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/indexing/Common.hpp"

#include "core/math/KahanSummation.h"
#include "core/DataTypes.h"

#ifdef DEBUG_ELEMENTWISE
#include "tinyhhg_core/format.hpp"
#endif

namespace hhg {
  namespace vertexdof {
    namespace macroedge {

      using walberla::real_c;
      using indexing::Index;

      inline void fillLocalCoords( const uint_t & level, uint_t i, const std::array< stencilDirection, 3>& element,
                                   const std::array<real_t*, 2>& coords, real_t localCoords[6] )
      {
        localCoords[0] = coords[0][vertexdof::macroedge::indexFromVertex( level, i, element[0] )];
        localCoords[1] = coords[1][vertexdof::macroedge::indexFromVertex( level, i, element[0] )];
        localCoords[2] = coords[0][vertexdof::macroedge::indexFromVertex( level, i, element[1] )];
        localCoords[3] = coords[1][vertexdof::macroedge::indexFromVertex( level, i, element[1] )];
        localCoords[4] = coords[0][vertexdof::macroedge::indexFromVertex( level, i, element[2] )];
        localCoords[5] = coords[1][vertexdof::macroedge::indexFromVertex( level, i, element[2] )];

#ifdef DEBUG_ELEMENTWISE
        WALBERLA_LOG_INFO_ON_ROOT( "Vertex Indices: " 
                                   << vertexdof::macroedge::indexFromVertex( level, i, element[0] ) << ", "
                                   << vertexdof::macroedge::indexFromVertex( level, i, element[1] ) << ", "
                                   << vertexdof::macroedge::indexFromVertex( level, i, element[2] ) );

        WALBERLA_LOG_INFO_ON_ROOT( "Local Coordinates: "
                                   << "(" << localCoords[0] << ", " << localCoords[1] << ") -- "
                                   << "(" << localCoords[2] << ", " << localCoords[3] << ") -- "
                                   << "(" << localCoords[4] << ", " << localCoords[5] << ")\n" );
#endif

      }


      inline void assembleStencilForMicroNode( const uint_t& level, uint_t i, Edge &edge, 
                                               const std::array<real_t*, 2>& globalCoords,
                                               std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                               std::vector<real_t>& edgeStencil )
      {

        using namespace P1Elements;

        typedef std::array< stencilDirection, 3 > Element;
        Element elementSW = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_S  }};
        Element elementS  = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_S,  stencilDirection::VERTEX_SE }};
        Element elementSE = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E  }};
        Element elementNE = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_E,  stencilDirection::VERTEX_N  }};
        Element elementN  = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_N,  stencilDirection::VERTEX_NW }};
        Element elementNW = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_W  }};

        real_t localCoords[6];
        Matrix3r localStiffness;

        std::fill(edgeStencil.begin(), edgeStencil.end(), walberla::real_c(0.0));

        fillLocalCoords( level, i, elementSW, globalCoords, localCoords);
        computeElementMatrix(localStiffness, localCoords);
        assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSW ), {{0,1,2}}, localStiffness, edgeStencil);

        fillLocalCoords( level, i, elementS, globalCoords, localCoords);
        computeElementMatrix(localStiffness, localCoords);
        assembleP1LocalStencil( convertStencilDirectionsToIndices( elementS ), {{0,1,2}}, localStiffness, edgeStencil);

        fillLocalCoords( level, i, elementSE, globalCoords, localCoords);
        computeElementMatrix(localStiffness, localCoords);
        assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSE ), {{0,1,2}}, localStiffness, edgeStencil);

        if (edge.getNumNeighborFaces() == 2)
          {
            fillLocalCoords( level, i, elementNE, globalCoords, localCoords);
            computeElementMatrix(localStiffness, localCoords);
            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNE ), {{0,1,2}}, localStiffness, edgeStencil);

            fillLocalCoords( level, i, elementN, globalCoords, localCoords);
            computeElementMatrix(localStiffness, localCoords);
            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementN ), {{0,1,2}}, localStiffness, edgeStencil);

            fillLocalCoords( level, i, elementNW, globalCoords, localCoords);
            computeElementMatrix(localStiffness, localCoords);
            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNW ), {{0,1,2}}, localStiffness, edgeStencil);
          }
      }


      template< typename ValueType >
      inline void applyElementwise( const uint_t & level, Edge &edge,
                                    const std::shared_ptr< PrimitiveStorage >& storage,
                                    std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                    std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Edge>, 2> &coordIds,
                                    UpdateType update ) {

        using namespace P1Elements;

        uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
        uint_t inner_rowsize = rowsize;

        auto src = edge.getData(srcId)->getPointer(level);
        auto dst = edge.getData(dstId)->getPointer(level);
        std::array<ValueType*, 2> globalCoords{{edge.getData(coordIds[0])->getPointer(level),
              edge.getData(coordIds[1])->getPointer(level)}};

        // typedef std::array< stencilDirection, 3 > Element;
        // Element elementSW = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_S  }};
        // Element elementS  = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_S,  stencilDirection::VERTEX_SE }};
        // Element elementSE = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E  }};
        // Element elementNE = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_E,  stencilDirection::VERTEX_N  }};
        // Element elementN  = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_N,  stencilDirection::VERTEX_NW }};
        // Element elementNW = {{ stencilDirection::VERTEX_C, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_W  }};
        // real_t localCoords[6];
        // Matrix3r localStiffness;

        ValueType tmp;
        std::vector<real_t> edgeStencil(7);
  
        for (size_t i = 1; i < rowsize - 1; ++i) {

          // std::fill(edgeStencil.begin(), edgeStencil.end(), walberla::real_c(0.0));

//           fillLocalCoords( level, i, elementSW, globalCoords, localCoords);
//           computeElementMatrix(localStiffness, localCoords);
//           assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSW ), {{0,1,2}}, localStiffness, edgeStencil);

//           fillLocalCoords( level, i, elementS, globalCoords, localCoords);
//           computeElementMatrix(localStiffness, localCoords);
//           assembleP1LocalStencil( convertStencilDirectionsToIndices( elementS ), {{0,1,2}}, localStiffness, edgeStencil);

//           fillLocalCoords( level, i, elementSE, globalCoords, localCoords);
//           computeElementMatrix(localStiffness, localCoords);
//           assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSE ), {{0,1,2}}, localStiffness, edgeStencil);

//           if (edge.getNumNeighborFaces() == 2)
//             {
//               fillLocalCoords( level, i, elementNE, globalCoords, localCoords);
//               computeElementMatrix(localStiffness, localCoords);
//               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNE ), {{0,1,2}}, localStiffness, edgeStencil);

//               fillLocalCoords( level, i, elementN, globalCoords, localCoords);
//               computeElementMatrix(localStiffness, localCoords);
//               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementN ), {{0,1,2}}, localStiffness, edgeStencil);

//               fillLocalCoords( level, i, elementNW, globalCoords, localCoords);
//               computeElementMatrix(localStiffness, localCoords);
//               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNW ), {{0,1,2}}, localStiffness, edgeStencil);
//             }

// #ifdef DEBUG_ELEMENTWISE
//           WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format("Edge.id = %d", edge.getID().getID() ));
//           for( uint_t weight = 0; weight < 7; ++weight )
//             {
//               WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, edgeStencil[weight] ) );
//             }
// #endif

          assembleStencilForMicroNode( level, i, edge, globalCoords, computeElementMatrix, edgeStencil );

          tmp = edgeStencil[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ]
            * src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

          for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
            {
              tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ]
                * src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
            }

          for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
            {
              tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ]
                * src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
            }

          if (edge.getNumNeighborFaces() == 2)
            {
              for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
                {
                  tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ]
                    * src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
                }
            }

          if (update == Replace) {
            dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
          } else if (update == Add) {
            dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += tmp;
          }
        }
      }

    } // namespace macroedge
  } // namespace vertexdof
} // namespace hhg

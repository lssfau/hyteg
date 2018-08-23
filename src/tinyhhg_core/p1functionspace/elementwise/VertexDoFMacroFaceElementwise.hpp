#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/indexing/Common.hpp"

#ifdef DEBUG_ELEMENTWISE
#include "tinyhhg_core/format.hpp"
#endif

namespace hhg {
  namespace vertexdof {
    namespace macroface {

      using walberla::uint_t;
      using walberla::real_c;
      using indexing::Index;


      inline void fillLocalCoords( const uint_t & Level, uint_t i, uint_t j, const P1Elements::P1Elements2D::P1Element& element,
                                   const std::array<real_t*, 2>& coords, real_t localCoords[6]) {

        localCoords[0] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[0] )];
        localCoords[1] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[0] )];
        localCoords[2] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[1] )];
        localCoords[3] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[1] )];
        localCoords[4] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[2] )];
        localCoords[5] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[2] )];
      }


      inline void assembleStencilForMicroNode( const uint_t& level, uint_t i, uint_t j, const std::array<real_t*, 2>& globalCoords,
                                               std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                               std::vector<real_t>& faceStencil )
      {

        using namespace P1Elements::P1Elements2D;

        real_t localCoords[6];
        Matrix3r localStiffness;

        std::fill( faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0) );
 
        for( uint_t k = 0; k < P1GrayElements.size(); ++k ) {

          // fill local coords
          localCoords[0] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][0] )];
          localCoords[1] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][0] )];
          localCoords[2] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][1] )];
          localCoords[3] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][1] )];
          localCoords[4] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][2] )];
          localCoords[5] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1GrayElements[k][2] )];

          // compute stencil contributions
          computeElementMatrix( localStiffness, localCoords );
          assembleP1LocalStencil( P1GrayStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil );
        }

        for (uint_t k = 0; k < P1BlueElements.size(); ++k) {

          // fill local coords
          localCoords[0] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][0] )];
          localCoords[1] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][0] )];
          localCoords[2] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][1] )];
          localCoords[3] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][1] )];
          localCoords[4] = globalCoords[0][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][2] )];
          localCoords[5] = globalCoords[1][vertexdof::macroface::indexFromVertex( level, i, j, P1BlueElements[k][2] )];

          // compute stencil contributions
          computeElementMatrix( localStiffness, localCoords );
          assembleP1LocalStencil( P1BlueStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil );
        }

#ifdef DEBUG_ELEMENTWISE
        WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format("FACE.id = %d", face.getID().getID() ));
        for( uint_t weight = 0; weight < 7; ++weight )
          {
            WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, faceStencil[weight] ) );
          }
#endif
      }


      template< typename ValueType >
      inline void applyElementwise( const uint_t & level, Face &face,
                                    std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                    std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Face>, 2> &coordIds,
                                    UpdateType update ) {

        using namespace P1Elements::P1Elements2D;
        typedef stencilDirection SD;

        uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
        uint_t inner_rowsize = rowsize;

        auto src = face.getData(srcId)->getPointer(level);
        auto dst = face.getData(dstId)->getPointer(level);
        std::array<ValueType*, 2> globalCoords{{face.getData(coordIds[0])->getPointer(level),
              face.getData(coordIds[1])->getPointer(level)}};

        ValueType tmp;
        // real_t localCoords[6];
        // Matrix3r localStiffness;
        std::vector<real_t> faceStencil(7);

        for (uint_t j = 1; j < rowsize - 2; ++j) {
          for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

//             std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));
 
//             for (uint_t k = 0; k < FaceVertexDoF::P1GrayElements.size(); ++k) {

//               // fill local coords
//               fillLocalCoords( level, i, j, FaceVertexDoF::P1GrayElements[k], globalCoords, localCoords);

//               // compute stencil
//               computeElementMatrix(localStiffness, localCoords);
//               assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
//             }

//             for (uint_t k = 0; k < FaceVertexDoF::P1BlueElements.size(); ++k) {

//               // fill local coords
//               fillLocalCoords( level, i, j, FaceVertexDoF::P1BlueElements[k], globalCoords, localCoords);

//               // fill coords
//               computeElementMatrix(localStiffness, localCoords);
//               assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
//             }

// #ifdef DEBUG_ELEMENTWISE
//             WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format("FACE.id = %d", face.getID().getID() ));
//             for( uint_t weight = 0; weight < 7; ++weight )
//               {
//                 WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, faceStencil[weight] ) );
//               }
// #endif

            assembleStencilForMicroNode( level, i, j, globalCoords, computeElementMatrix, faceStencil );

            if (update == Replace) {
              tmp = ValueType(0);
            }
            else {
              tmp = dst[vertexdof::macroface::indexFromVertex( level, i, j, SD::VERTEX_C )];
            }

            tmp += faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, SD::VERTEX_C )];

            // strangely the intel compiler can't handle this if it is a loop
            static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );

            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]
              * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];

            dst[vertexdof::macroface::indexFromVertex( level, i, j, SD::VERTEX_C )] = tmp;
          }
          --inner_rowsize;
        }
      }


      template< typename ValueType >
      inline void point_smooth_elementwise( const uint_t &level, Face &face,
                                            std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                            std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Face>, 2> &coordIds,
                                            ValueType relax ) {

        uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
        uint_t inner_rowsize = rowsize;

        auto dst = face.getData(dstId)->getPointer( level );
        auto rhs = face.getData(rhsId)->getPointer( level );
        auto src = face.getData(srcId)->getPointer( level );

        ValueType tmp;
        std::vector<real_t> faceStencil(7);
        std::array<ValueType*, 2> globalCoords{{face.getData(coordIds[0])->getPointer(level),
              face.getData(coordIds[1])->getPointer(level)}};

        for( uint_t j = 1; j < rowsize - 2; ++j ) {
          for( uint_t i = 1; i < inner_rowsize - 2; ++i ) {

            assembleStencilForMicroNode( level, i, j, globalCoords, computeElementMatrix, faceStencil );

            tmp = rhs[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )];

            for( uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k ) {

              tmp -= faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]
                * src[vertexdof::macroface::indexFromVertex( level, i, j, vertexdof::macroface::neighborsWithoutCenter[k] )];
            }

            dst[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] =
              ( ValueType(1.0) - relax ) * src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )]
              + relax * tmp / faceStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)];
          }
          --inner_rowsize;
        }
      }

    } // namespace macroface
  } // namespace vertexdof
} // namespace hhg

#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
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


      inline void fillLocalCoords( const uint_t & Level, uint_t i, uint_t j, const P1Elements::P1Element& element,
                                   const std::array<real_t*, 2>& coords, real_t localCoords[6]) {

        localCoords[0] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[0] )];
        localCoords[1] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[0] )];
        localCoords[2] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[1] )];
        localCoords[3] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[1] )];
        localCoords[4] = coords[0][vertexdof::macroface::indexFromVertex( Level, i, j, element[2] )];
        localCoords[5] = coords[1][vertexdof::macroface::indexFromVertex( Level, i, j, element[2] )];
      }


      template< typename ValueType >
      inline void applyElementwise( const uint_t & Level, Face &face,
                                    std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                    std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Face>, 2> &coordIds,
                                    UpdateType update ) {

        using namespace P1Elements;
        typedef stencilDirection SD;

        uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
        uint_t inner_rowsize = rowsize;

        auto src = face.getData(srcId)->getPointer(Level);
        auto dst = face.getData(dstId)->getPointer(Level);
        std::array<ValueType*, 2> globalCoords{{face.getData(coordIds[0])->getPointer(Level),
              face.getData(coordIds[1])->getPointer(Level)}};

        ValueType tmp;
        real_t localCoords[6];
        Matrix3r localStiffness;
        std::vector<real_t> faceStencil(7);

        for (uint_t j = 1; j < rowsize - 2; ++j) {
          for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

            std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));
 
            for (uint_t k = 0; k < FaceVertexDoF::P1GrayElements.size(); ++k) {

              // fill local coords
              fillLocalCoords( Level, i, j, FaceVertexDoF::P1GrayElements[k], globalCoords, localCoords);

              // compute stencil
              computeElementMatrix(localStiffness, localCoords);
              assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
            }

            for (uint_t k = 0; k < FaceVertexDoF::P1BlueElements.size(); ++k) {

              // fill local coords
              fillLocalCoords( Level, i, j, FaceVertexDoF::P1BlueElements[k], globalCoords, localCoords);

              // fill coords
              computeElementMatrix(localStiffness, localCoords);
              assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
            }

#ifdef DEBUG_ELEMENTWISE
            WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format("FACE.id = %d", face.getID().getID() ));
            for( uint_t weight = 0; weight < 7; ++weight )
              {
                WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, faceStencil[weight] ) );
              }
#endif

            if (update == Replace) {
              tmp = ValueType(0);
            }
            else {
              tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
            }

            tmp += faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];

            // strangely the intel compiler can't handle this if it is a loop
            static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );

            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
            tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];

            dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;
          }
          --inner_rowsize;
        }
      }

    } // namespace macroface
  } // namespace vertexdof
} // namespace hhg

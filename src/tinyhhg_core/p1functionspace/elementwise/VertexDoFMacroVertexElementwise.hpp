#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

#ifdef DEBUG_ELEMENTWISE
#include "tinyhhg_core/format.hpp"
#endif

namespace hhg {
  namespace vertexdof {
    namespace macrovertex {

      inline void fillLocalCoords( const std::array<uint_t, 3>& element, const std::array<real_t*, 2>& coords,
                                   real_t localCoords[6] ) {

        localCoords[0] = coords[0][element[0]];
        localCoords[1] = coords[1][element[0]];
        localCoords[2] = coords[0][element[1]];
        localCoords[3] = coords[1][element[1]];
        localCoords[4] = coords[0][element[2]];
        localCoords[5] = coords[1][element[2]];
      }


      inline void assembleStencilForMicroNode( const uint_t& level, Vertex &vertex, 
                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                               const std::array<real_t*, 2>& globalCoords,
                                               std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                               std::vector<real_t>& vertexStencil ) {

        using namespace P1Elements;

        Matrix3r localStiffness;
        real_t localCoords[6];

        // iterate over adjacent faces
        for( auto& faceId : vertex.neighborFaces() ) {

          Face* face = storage->getFace(faceId);

          uint_t v_i = face->vertex_index(vertex.getID());

          std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

          std::array<uint_t, 3> stencilMap;
          stencilMap[0] = 0;

          std::array<uint_t, 3> dofMap;
          dofMap[0] = v_i;

          // iterate over adjacent edges
          for (uint_t i = 0; i < adj_edges.size(); ++i) {
            uint_t edge_idx = vertex.edge_index(adj_edges[i]) + 1;
            Edge* edge = storage->getEdge(adj_edges[i]);
            PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

            stencilMap[i+1] = edge_idx;
            dofMap[i+1] = face->vertex_index(vertex_j);
          }

          fillLocalCoords(stencilMap, globalCoords, localCoords);
          computeElementMatrix(localStiffness, localCoords);

          assembleP1LocalStencil(stencilMap, {{0,1,2}}, localStiffness, vertexStencil);
        }

#ifdef DEBUG_ELEMENTWISE
        WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( "Vertex.id = %d", vertex.getID().getID() ));
        for( uint_t weight = 0; weight < vertexStencil.size(); ++weight )
          {
            WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, vertexStencil[weight] ) );
          }
#endif
      }


      template< typename ValueType>
      inline void applyElementwise( uint_t level, Vertex &vertex,
                                    const std::shared_ptr< PrimitiveStorage >& storage,
                                    std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                                    std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Vertex>, 2> &coordIds,
                                    UpdateType update ) {

        using namespace P1Elements;

        auto src = vertex.getData(srcId)->getPointer(level);
        auto dst = vertex.getData(dstId)->getPointer(level);
        std::array<ValueType*, 2> globalCoords{{vertex.getData(coordIds[0])->getPointer(level),
              vertex.getData(coordIds[1])->getPointer(level)}};

        std::vector<real_t> vertexStencil(1 + vertex.getNumNeighborEdges());

//         Matrix3r localStiffness;
//         real_t localCoords[6];

//         // iterate over adjacent faces
//         for (auto& faceId : vertex.neighborFaces()) {

//           Face* face = storage->getFace(faceId);

//           uint_t v_i = face->vertex_index(vertex.getID());

//           std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

//           std::array<uint_t, 3> stencilMap;
//           stencilMap[0] = 0;

//           std::array<uint_t, 3> dofMap;
//           dofMap[0] = v_i;

//           // iterate over adjacent edges
//           for (uint_t i = 0; i < adj_edges.size(); ++i) {
//             uint_t edge_idx = vertex.edge_index(adj_edges[i]) + 1;
//             Edge* edge = storage->getEdge(adj_edges[i]);
//             PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

//             stencilMap[i+1] = edge_idx;
//             dofMap[i+1] = face->vertex_index(vertex_j);
//           }

//           fillLocalCoords(stencilMap, globalCoords, localCoords);
//           computeElementMatrix(localStiffness, localCoords);

//           assembleP1LocalStencil(stencilMap, {{0,1,2}}, localStiffness, vertexStencil);
//         }

// #ifdef DEBUG_ELEMENTWISE
//         WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( "Vertex.id = %d", vertex.getID().getID() ));
//         for( uint_t weight = 0; weight < vertexStencil.size(); ++weight )
//           {
//             WALBERLA_LOG_DEVEL_ON_ROOT( hhg::format( " Stencil weight[%d] = %e", weight, vertexStencil[weight] ) );
//           }
// #endif

        assembleStencilForMicroNode( level, vertex, storage, globalCoords, computeElementMatrix, vertexStencil );

        if (update==Replace) {
          dst[0] = vertexStencil[0]*src[0];
        } else if (update==Add) {
          dst[0] += vertexStencil[0]*src[0];
        }

        for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
          dst[0] += vertexStencil[i + 1]*src[i + 1];
        }
      }

    } // namespace macrovertex
  } // namespace vertexdof
} // namespace hhg

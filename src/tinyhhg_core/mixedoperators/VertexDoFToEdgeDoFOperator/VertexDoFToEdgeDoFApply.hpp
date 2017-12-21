#pragma once

#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg{
namespace VertexDoFToEdgeDoF{


template<uint_t Level>
inline void applyFaceTmpl(Face &face,
                          const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                          UpdateType update){

}

SPECIALIZE(void, applyFaceTmpl, applyFace)

template<uint_t Level>
inline void applyEdgeTmpl(Edge &edge,
                          const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                          UpdateType update) {
}


SPECIALIZE(void, applyEdgeTmpl, applyEdge)

} // VertexDoFToEdgeDoFOperator
} // namespace hhg
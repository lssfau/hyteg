#pragma once

#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "VertexDoFToEdgeDoFMemory.hpp"



namespace hhg{
namespace VertexDoFToEdgeDoF{

template<size_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceVertexDoFToEdgeDoFStencilMemory< real_t >, Face> &operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId, UpdateType update) {

}

} // VertexDoFToEdgeDoFOperator
} // namespace hhg
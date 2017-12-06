#pragma once

#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg{
namespace EdgeDoFToVertexDoFVertex{

inline void apply(uint_t level,
                  Vertex &vertex,
                  const PrimitiveDataID<StencilMemory < real_t >, Vertex> &operatorId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &srcId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstId,
                  UpdateType update)
{
  WALBERLA_LOG_DEVEL("TODO")
}

}

namespace EdgeDoFToVertexDoFEdge{

template<uint_t Level>
inline void apply_tmpl(Edge &edge,
                  const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                  UpdateType update)
{
  WALBERLA_LOG_DEVEL("TODO")
}

SPECIALIZE(void, apply_tmpl, apply)

}

namespace EdgeDoFToVertexDoFFace{

template<uint_t Level>
inline void apply_tmpl(Face &face,
                       const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                       UpdateType update)
{
  WALBERLA_LOG_DEVEL("TODO")
}

SPECIALIZE(void, apply_tmpl, apply)

}


}
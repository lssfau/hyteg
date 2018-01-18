#pragma once

#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"


namespace hhg{
namespace P2{

namespace vertex {

void smoothGSvertexDof(Vertex vertex,
                       const PrimitiveDataID<StencilMemory< real_t >, Vertex> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Vertex> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &getVertexDoFID,
                       uint_t level);

} /// namespace vertex

namespace edge {

void smoothGSvertexDof(Edge edge,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsVertexDoFID,
                       uint_t level);

void smoothGSedgeDof(Edge edge,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsEdgeDoFID,
                       uint_t level);

} /// namespace edge

namespace face {

void smoothGSvertexDof(Face face,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &vertexDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &edgeDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &rhsVertexDoFID,
                     uint_t level);

void smoothGSedgeDof(Face face,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &vertexDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &edgeDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &rhsEdgeDoFID,
                     uint_t level);

} /// namespace face

} /// namespace P2
} /// namespace hhg
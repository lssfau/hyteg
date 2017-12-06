#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {

class EdgeDoFToVertexDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::P1Function < real_t > >
{
public:

  EdgeDoFToVertexDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel);
  ~EdgeDoFToVertexDoFOperator() final = default;

  void apply_impl(EdgeDoFFunction< real_t >& src,  P1Function< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID<StencilMemory< real_t >, Vertex> &getVertexStencilID_() const;
  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID_() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID_() const;

private:
  PrimitiveDataID<StencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;

};

}// namespace hhg
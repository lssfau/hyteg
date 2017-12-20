#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include "EdgeDoFFunction.hpp"

namespace hhg{

class EdgeDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::EdgeDoFFunction < real_t > >
{
public:

  EdgeDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel);
  ~EdgeDoFOperator() final = default;

  void apply_impl(EdgeDoFFunction< real_t >& src,  EdgeDoFFunction< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID_() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID_() const;

private:
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;

};


}
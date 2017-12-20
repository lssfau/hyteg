#include "EdgeDoFOperator.hpp"

namespace hhg {

EdgeDoFOperator::EdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                 size_t minLevel,
                                 size_t maxLevel)
  : Operator(storage, minLevel, maxLevel)
{
  WALBERLA_ABORT("implement me");
}

void
EdgeDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src, EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) {
  WALBERLA_ABORT("implement me");
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID_() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID_() const {
  return faceStencilID_;
}

}


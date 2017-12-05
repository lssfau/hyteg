#include "EdgeDoFToVertexDoFOperator.hpp"

namespace hhg{

EdgeDoFToVertexDoFOperator::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                                       size_t minLevel,
                                                       size_t maxLevel)
  :Operator(storage,minLevel,maxLevel)
{
  WALBERLA_ABORT("implement me");
}

void EdgeDoFToVertexDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src,
                                            P1Function<real_t> &dst,
                                            uint_t level,
                                            DoFType flag,
                                            UpdateType updateType)
{
  WALBERLA_ABORT("implement me");
}

const PrimitiveDataID<VertexEdgeDoFToVertexDoFStencilMemory < real_t >, Vertex > &EdgeDoFToVertexDoFOperator::getVertexStencilID_() const {
  return vertexStencilID_;
}

const PrimitiveDataID<real_t> &EdgeDoFToVertexDoFOperator::getEdgeStencilID_() const {
  return edgeStencilID_;
}

const PrimitiveDataID<real_t> &EdgeDoFToVertexDoFOperator::getFaceStencilID_() const {
  return faceStencilID_;
}


}
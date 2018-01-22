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

  real_t* getFaceStencil(const PrimitiveID& faceId, uint_t level);

  void apply_impl(EdgeDoFFunction< real_t >& src,  EdgeDoFFunction< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID() const;

private:
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;

};

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);


}
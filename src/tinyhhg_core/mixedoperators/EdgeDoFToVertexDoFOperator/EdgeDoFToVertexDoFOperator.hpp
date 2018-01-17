#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {

class EdgeDoFToVertexDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::P1Function < real_t > >
{
public:

  EdgeDoFToVertexDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, const size_t & minLevel, const size_t & maxLevel);
  ~EdgeDoFToVertexDoFOperator() final = default;

  real_t* getVertexStencil(const PrimitiveID& vertexId, uint_t level);
  real_t* getEdgeStencil(const PrimitiveID& edgeId, uint_t level);
  real_t* getFaceStencil(const PrimitiveID& faceId, uint_t level);

  void apply_impl(EdgeDoFFunction< real_t >& src,  P1Function< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID<StencilMemory< real_t >, Vertex> &getVertexStencilID() const;
  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID() const;

private:
  PrimitiveDataID<StencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;

};

namespace EdgeDoFToVertexDoF {

/// \param level the stencil size is independent of the level
/// \param numDependencies number of adjacent edges of the vertex
/// \return number of entries in the stencil on a macro vertex, be aware that this one entry to large on the boundaries
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// \param level the stencil size is independent of the level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of entries in the stencil on a macro edge
uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// \param level the stencil size is independent of the level
/// \param numDependencies on a macro face this is always constant
/// \return number of entries in the stencil on a macro face
uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

}/// namespace EdgeDoFToVertexDoF
}/// namespace hhg
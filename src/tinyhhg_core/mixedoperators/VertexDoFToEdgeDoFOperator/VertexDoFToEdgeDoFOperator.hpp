#pragma once

#include "VertexDoFToEdgeDoFApply.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif


namespace hhg {

class VertexDoFToEdgeDoFOperator : public Operator<P1Function< real_t >, EdgeDoFFunction< real_t > >
{
public:
  VertexDoFToEdgeDoFOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel);
  ~VertexDoFToEdgeDoFOperator() final = default;

  void apply_impl(P1Function< real_t > & src, EdgeDoFFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace);

  /// since the Vertex does not own any EdgeDoFs only edge and face are needed
  const PrimitiveDataID< StencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }
  const PrimitiveDataID< StencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }

private:
  PrimitiveDataID< StencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID< StencilMemory< real_t >, Face> faceStencilID_;

};

namespace VertexDoFToEdgeDoF {

/// \param level stencil size is independent of level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// \param level stencil size is independent of level
/// \param numDependencies not needed for faces
/// \return number of the stencil entries
uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);
}

}/// namespace hhg

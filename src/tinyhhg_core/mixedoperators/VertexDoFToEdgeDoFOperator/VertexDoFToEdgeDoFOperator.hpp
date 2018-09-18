#pragma once

#include "VertexDoFToEdgeDoFApply.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p2functionspace/generated/p2_divt.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif


namespace hhg {

template<class UFCOperator>
class VertexDoFToEdgeDoFOperator : public Operator<P1Function< real_t >, EdgeDoFFunction< real_t > >
{
public:
  VertexDoFToEdgeDoFOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel);
  ~VertexDoFToEdgeDoFOperator() final = default;

  void apply_impl(P1Function< real_t > & src, EdgeDoFFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace);

  /// since the Vertex does not own any EdgeDoFs only edge, face and cell are needed
  const PrimitiveDataID< StencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }
  const PrimitiveDataID< StencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }
  const PrimitiveDataID< StencilMemory< real_t >, Cell> &getCellStencilID() const { return cellStencilID_; }

private:
  void assembleStencils();
  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type);

  PrimitiveDataID< StencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID< StencilMemory< real_t >, Face> faceStencilID_;
  PrimitiveDataID< StencilMemory< real_t >, Cell> cellStencilID_;

};

namespace VertexDoFToEdgeDoF {

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroCellVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );
}

typedef VertexDoFToEdgeDoFOperator<hhg::fenics::NoAssemble> GenericVertexDoFToEdgeDoFOperator;
typedef VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_0_otherwise> VertexToEdgeDivTxOperator;
typedef VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_1_otherwise> VertexToEdgeDivTyOperator;

}/// namespace hhg

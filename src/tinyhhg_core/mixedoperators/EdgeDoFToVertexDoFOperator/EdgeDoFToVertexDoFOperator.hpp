#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p2functionspace/generated/p2_div.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hhg {

template<class UFCOperator>
class EdgeDoFToVertexDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::P1Function < real_t > >
{
public:

  EdgeDoFToVertexDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, const size_t & minLevel, const size_t & maxLevel);
  ~EdgeDoFToVertexDoFOperator() final = default;

  void apply_impl(EdgeDoFFunction< real_t >& src,  P1Function< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID< StencilMemory< real_t >, Vertex> &getVertexStencilID() const;
  const PrimitiveDataID< StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID< StencilMemory< real_t >, Face  > &getFaceStencilID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face  > &getFaceStencil3DID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell  > &getCellStencilID() const;

private:
  void assembleStencils();
  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type);

  PrimitiveDataID<StencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face  > faceStencil3DID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell  > cellStencilID_;

};

namespace EdgeDoFToVertexDoF {

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro vertex, be aware that this one entry to large on the boundaries
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro edge
uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro face
uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro cell
uint_t macroCellEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

}/// namespace EdgeDoFToVertexDoF

typedef EdgeDoFToVertexDoFOperator<hhg::fenics::NoAssemble> GenericEdgeDoFToVertexDoFOperator;
typedef EdgeDoFToVertexDoFOperator<p2_div_cell_integral_0_otherwise> EdgeToVertexDivxOperator;
typedef EdgeDoFToVertexDoFOperator<p2_div_cell_integral_1_otherwise> EdgeToVertexDivyOperator;

}/// namespace hhg

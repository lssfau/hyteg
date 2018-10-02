#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroCell.hpp"

#include "EdgeDoFFunction.hpp"

namespace hhg{

class EdgeDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::EdgeDoFFunction < real_t > >
{
public:

  EdgeDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel);
  ~EdgeDoFOperator() final = default;

  void apply_impl(EdgeDoFFunction< real_t >& src,  EdgeDoFFunction< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID() const;
  const PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell  > &getCellStencilID() const;

private:
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;
  PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell  > cellStencilID_;

};

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

uint_t macroCellEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

}
